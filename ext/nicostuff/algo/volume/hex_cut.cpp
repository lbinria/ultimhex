#include <algo/volume/hex_cut.h>

#include <ultimaille/all.h>
#include <framework/trace.h>
#include "toolbox.h"
#include "drop_attribute.h"

#include <algo/volume/hex_select.h>
#include <algo/volume/hex_edit.h>
#include <algo/volume/hex_smooth.h>



HexCutter::HexCutter(Hexahedra& hex) : hex(hex), chart(hex, -1){
	//um_assert(hex.connected());
	drop_all_steps=true;
}

HexCutter::~HexCutter() {
	for (auto s : sdf) delete s;
}

void HexCutter::init_from_charts() {
	int nb_charts = 0;
	FOR(f, hex.nfacets())  nb_charts = std::max(nb_charts, chart[f] + 1);
	FOR(c, nb_charts) sdf.push_back(new TrianglesSDF(hex, chart, c));

}

void HexCutter::optimize_SDF(int sdf_id) {
	
	HexPointSelect variable_vertices(hex);
	variable_vertices.fill(false);
	variable_vertices.by_facet_value(chart, sdf_id, true);
	FOR(iter, 3) variable_vertices.dilate(true);
	//variable_vertices.drop();

	HexSmoother smoother(hex);
	smoother.set_lock(variable_vertices, false);


	Drop(hex,chart).apply("charts");
	for (auto h : hex.iter_halfedges()) if (chart[h.facet()] != -1) 
//if (chart[h.facet()] == 	0//sdf_id) 
{
		auto sdf_ptr = sdf[chart[h.facet()]];
		vec3 pos = h.from().pos();
		vec3 P = sdf_ptr->proj(pos);
		vec3 n = sdf_ptr->grad(P);
		smoother.add_constraint(h.from(), P, n);
	}


	////// begin test sdf proj
	PolyLine pl;
	pl.create_edges(smoother.tan_constraints.size());
	pl.points.create_points(2*smoother.tan_constraints.size());
	FOR(i,smoother.tan_constraints.size()){
		FOR(lv,2) pl.vert(i,lv)=2*i+lv;
		vec3 v0 = hex.points[smoother.tan_constraints[i].v];
		vec3 v1 = smoother.tan_constraints[i].P;
		pl.points[2*i] = v0;
		pl.points[2*i+1] = v1;
	}
	DropPolyLine(pl).apply("proj");
	DropSurface(dynamic_cast<TrianglesSDF*> (sdf[0])->m).apply("sdf");
	/////  end test sdf proj

	Trace::step("LS opti");
	//smoother.show_scaled_jacobien("SJ before LS", 0);
	 //Drop(hex, chart)._select([&](int f) { return chart[f] == sdf_id; }).apply(" OOOOOOOOOOOO New SDF LS OOOOOOOOO");
	smoother.smooth_LS();
	 //Drop(hex, chart)._select([&](int f) { return chart[f] == sdf_id; }).apply(" OOOOOOOOOOOO New SDF LS SMOOTHED OOOOOOOOO");
	//smoother.show_scaled_jacobien("SJ after LS", 0);

	//throw(std::exception("break"));


	if (drop_all_steps) Drop(hex, chart)._select([&](int f) { return chart[f] == sdf_id; }).apply(" OOOOOOOOOOOO New SDF LS smoothed OOOOOOOOO");
	if (drop_all_steps) variable_vertices.drop("variable vertices");
	if (drop_all_steps) smoother.show_scaled_jacobien("SJ after LS", 0);


	Trace::step("untangle");
	HexPointSelect badSJ(hex, false);
	CellAttribute<double> SJ(hex);
	ToolBox(hex).eval_quality(SJ);
	double minSJ = 1e20;
	for (auto c : hex.iter_cells()) minSJ = std::min(minSJ, SJ[c]);
	for (auto c : hex.iter_cells()) if (SJ[c] < 1.1 * std::max(0., minSJ)) badSJ.set_cell(c, true);
	if (minSJ < 0) {
		if (drop_all_steps) {
			badSJ.drop_cells_facet(chart, "constraints close to inverted cells");
			badSJ.drop_cells("inverted SJ", 8);
			badSJ.drop_cells("inverted SJ neig", 4);
		}
		FOR(l, 3) badSJ.dilate(true);
		//badSJ.drop("badSJ");

		//variable_vertices.drop("variable_vertices");

		for (auto v : hex.iter_vertices()) badSJ[v] = badSJ[v] && variable_vertices[v];
		badSJ.drop("untangle zone");
		//badSJ.invert();
		smoother.set_lock(badSJ,false);
		//if (drop_all_steps) smoother.show_scaled_jacobien("avt elliptic", 0);
		smoother.smooth_elliptic();
		if (drop_all_steps) {
			badSJ.drop_cells_facet(chart, "untangled constraints close to inverted cells");
			badSJ.drop_cells("untangled inverted SJ", 8);
			badSJ.drop_cells("untangled inverted SJ neig", 4);
		}

		if (drop_all_steps) smoother.show_scaled_jacobien("untangle failure", 0);
		//Drop(hex, chart)._shrink(0).apply("untangled");
	}

	Trace::step("Smooth");
	for (auto v : hex.iter_vertices()) badSJ[v] = !badSJ[v] && variable_vertices[v];
	smoother.set_lock(badSJ);
	smoother.smooth_elliptic();
	
	if (drop_all_steps) {
		variable_vertices.drop_cells_facet(chart, "constraints close to inverted cells");
		variable_vertices.drop_cells("modified hexes", 8);
	}

	//smoother.show_scaled_jacobien("ap smoothing", 1.1);
	//Drop(hex, chart)._shrink(0).apply("smoothed");
}




void HexCutter::show() {
	for (auto s : sdf) s->show(ToolBox(hex.points).bbox());
}

void HexCutter::show_feature_neigs(std::string name, bool include_boundary, bool only_last_SDF) {
	Hexahedra select;
	ToolBox(select.points).copy_from(hex.points);
	//EC3d conn(hex);
	int last_SDF = -1;
	for (auto f : hex.iter_facets())
		last_SDF = std::max(last_SDF, chart[f]);

	CellAttribute<bool> lock(hex, true);
	if (only_last_SDF) {
		for (auto f : hex.iter_facets())
			if (chart[f] == last_SDF) lock[f.cell()] = false;
	}
	else
		for (auto h : hex.iter_halfedges())
			if (chart[h.facet()] != -1) lock[h.cell()] = false;

	if (!include_boundary)
		for (auto h : hex.iter_halfedges())
			if (!h.opposite_c().active()) lock[h.cell()] = true;


	for (auto c : hex.iter_cells()) if (!lock[c]) {
		int nc = select.create_cells(1);
		FOR(lv, 8) select.vert(nc, lv) = hex.vert(c, lv);
	}
	select.delete_isolated_vertices();


	HexSmoother(select).show_scaled_jacobien(name);
}



void HexCutter::add_sdf(SDF* s) {
	int cur_chart_id = sdf.size();

	{// flag facets "on the iso-0" 

		// group blocks of same sign
		DisjointSet ds(hex.ncells());
		for (auto f : hex.iter_facets()) if (chart[f] != -1) {
			auto opp = f.opposite();
			if (!opp.active()) continue;
			ds.merge(f.cell(), opp.cell());
		}
		std::vector<int> cell2block;
		int nb_blocks = ds.get_sets_id(cell2block);

		// compute average signed distance of each block
		std::vector<double> block_dist(nb_blocks);
		for (auto c : hex.iter_cells())
			block_dist[cell2block[c]] += s->dist(ToolBox(hex).cell_geom(c).bary_verts());


		// make the iso_0 manifold 
		bool done = false;
		while (!done) {
			done = true;
			for (auto h : hex.iter_halfedges()) {
				bool border = false;
				int nb_switch = 0;
				for (auto cir : h.iter_CCW_around_edge()) {
					if (!cir.opposite_c().active()) border = true;
					else if (block_dist[cell2block[cir.cell()]] * block_dist[cell2block[cir.opposite_c().cell()]] < 0)
						nb_switch++;
					if (nb_switch > 2) {
						done = false;
						for (auto inner_cir : h.iter_CCW_around_edge())
							block_dist[cell2block[inner_cir.cell()]] = 1;
					}
				}
			}
		}



		// set chart that separates positive and negative SDF blocks
		for (auto f : hex.iter_facets()) {
			auto opp = f.opposite();
			if (!opp.active()) continue;
			if (block_dist[cell2block[f.cell()]] > 0 && block_dist[cell2block[opp.cell()]] <= 0) {
				chart[f] = cur_chart_id;
				chart[opp] = cur_chart_id;
			}
		}
		sdf.push_back(s);
	}

	{// pad the iso-0 facets

		CellFacetAttribute<bool> to_pad(hex, false);
		FOR(f, hex.nfacets())
			to_pad[f] = (chart[f] == cur_chart_id);
		
		Drop(hex,to_pad)._skip_value(false).apply("2pad");
		HexPad padder(hex);
		padder.apply(to_pad);

		hex.connect();
		for (auto f : hex.iter_facets()) if (chart[f] == cur_chart_id) {
			for (auto h : f.iter_halfedges()) {
				int old_chart = chart[h.opposite_f().facet()];
				chart[h.opposite_c().opposite_f().facet()] = old_chart;
			}
			chart[f] = -1;
		}
		for (auto f : hex.iter_facets()) if (to_pad[f])
			chart[f.halfedge(0).opposite_c().opposite_f().next().next().opposite_f().facet()] = cur_chart_id;
	}
	Drop(hex,chart).apply("hex sdf iso");
	optimize_SDF(cur_chart_id);

}

