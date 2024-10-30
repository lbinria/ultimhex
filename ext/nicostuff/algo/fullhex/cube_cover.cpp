#include "cube_cover.h"

#include <fullhex/gp_basic.h>

#include <framework/trace.h>

#include <drop_attribute.h>
#include <volume/frame3D.h>





using namespace UM::Linear;

struct CCParam {
	bool debug_GP_across_cut = true;
	bool debug_continuity_outside_cuts = true;
	bool debug_seamless_across_cuts = true;
	bool GP_on_singularities = true;
	bool align_with_boundary = true;
	bool snapped_to_boundary = true;
	bool enforce_locked_cell = true;
} paramCC;




CubeCover::CubeCover(Tetrahedra& p_m, FF3D& p_ff, CellCornerAttribute<vec3>& U, CellAttribute<bool>& locked_cells,CellFacetAttribute<int>& constraint_type) : m(p_m), ff(p_ff), U(U), locked_cells(locked_cells),constraint_type(constraint_type) {
um_assert(!"CubeCover is deprecarted");
}

void CubeCover::cut_graph_3D_and_brush_ff(CellFacetAttribute<bool>& cut) {
	Trace::Section sec("cut_graph_3D");
	int seed = -1;
	Trace::step("compute dist 2 feature");
	CellAttribute<double> dist(m, 1e20);
	{
		auto cmp = [&](int left, int right) { return dist[left] > dist[right]; };
		std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);

		for (auto f : m.iter_facets()) if (!f.opposite().active() && !locked_cells[f.cell()]) {
			if (std::abs(uvw_tri(U, f).normal()[2]) > .9) continue;
			dist[f.cell()] = 0;
			queue.push(f.cell());
		}

		while (!queue.empty()) {
			Volume::Cell c(m, queue.top());
			queue.pop();
			seed = c;
			for (auto f : c.iter_facets()) {
				if (!f.opposite().active()) continue;
				int opp = f.opposite().cell();
				if (dist[opp] > dist[c] + 1) {
					dist[opp] = dist[c] + 1;
					queue.push(opp);
				}
			}
		}
		Drop(m, dist).apply("dist");
	}
	for(auto c:m.iter_cells()) dist[c] = Tetrahedron(c).bary_verts()[0];
		Drop(m, dist).apply("dist");

	Trace::step("covering tree");
	{
		for (auto f : m.iter_facets()) cut[f] = true;
		CellAttribute<bool> visited(m, false);
		auto cmp = [&](int left, int right) { return dist[left] < dist[right]; };
		std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);
		queue.push(seed);
		visited[seed] = true;
		while (!queue.empty()) {
			Volume::Cell c(m, queue.top());
			queue.pop();
			for (auto f : c.iter_facets()) {
				auto opp_f = f.opposite();
				if (!opp_f.active()) continue;
				auto opp_c = opp_f.cell();
				if (visited[opp_c]) continue;
				AxisPermutation perm;
				perm = permute_Jj_to_approx_Ji(ff[c], ff[opp_c]);
				//FOR(lv, 4) U[opp_c.corner(lv)] = perm.get_mat() * U[opp_c.corner(lv)];

				ff[opp_c] = perm.get_mat() * ff[opp_c];

				cut[opp_f] = false;
				cut[f] = false;
				visited[opp_c] = true;
				queue.push(opp_c);
			}
		}
	}

	for (auto c:m.iter_cells()) if (locked_cells[c]){
		mat3x3 Ju=uvw_to_jacobian(c, U);
//		auto perm = permute_Jj_to_approx_Ji( Ju,ff[c]);
		auto perm = permute_Jj_to_approx_Ji( ff[c],Ju);
		//if (perm.mid!=0) plop(perm.mid);
		FOR(lv, 4) U[c.corner(lv)] = perm.get_mat() * U[c.corner(lv)];
	}

	ff.update_axis_permutation();

	Trace::step("remove boundaries and contractible loops ");

	for (auto f : m.iter_facets()) if (!f.opposite().active()) cut[f] = false;
	
	EdgeGraph eg(m);
	EdgeAttribute<int> nb_cuts(eg,0);
	EdgeAttribute<bool> on_border(eg,false);
	for(auto h:m.iter_halfedges()) if (cut[h.facet()])					nb_cuts [eg.edge_from_halfedge(h)]++;
	for(auto h:m.iter_halfedges()) if (!h.facet().opposite().active())	on_border[eg.edge_from_halfedge(h)]=true;
	bool done = false;
	while (!done) {
		done = true;
		for (auto f : m.iter_facets()) if (f.opposite().active() ) for(auto h:f.iter_halfedges()) if( cut[f]){
			if (nb_cuts[eg.edge_from_halfedge(h)]!=1) continue; 
			if (on_border[eg.edge_from_halfedge(h)]) continue; 
			if (ff.edge_is_singular(h)) continue;
			cut[f] = false;
			cut[f.opposite()] = false;
			for(auto cir : f.iter_halfedges())				nb_cuts[eg.edge_from_halfedge(cir)]--;
			for(auto cir : f.opposite().iter_halfedges())	nb_cuts[eg.edge_from_halfedge(cir)]--;
			done = false;
		}
	}

	{
		Trace::step("render cut graph");
		TetBoundary tb(m);
		DropSurface(tb.tri).apply("border");
		Drop(m, cut)._skip_value(false).apply("cut_graph");
		EdgeAttribute<int> cut_graph_border(eg,0);
		for(auto h:m.iter_halfedges()) if (!h.facet().opposite().active()) cut_graph_border[eg.edge_from_halfedge(h)]=nb_cuts[eg.edge_from_halfedge(h)];
		Drop<PolyLine,EdgeAttribute<int> >(eg,cut_graph_border)._skip_value(0).apply_wireframe("cut_graph_border");
	}

}

void CubeCover::add_real_constraints(CellFacetAttribute<bool>& cut, ConstrainedLeastSquares& ls) {
	
	Trace::Section sec("real constraints");
	
	if (paramCC.debug_continuity_outside_cuts) {
		Trace::step("debug_continuity_outside_cuts");
		for (auto f : m.iter_facets()) if (!cut[f]) FOR(lh, 3) {
			auto h = f.halfedge(lh);
			auto opp_h = h.opposite_c();
			if (!opp_h.active()) continue;
			FOR(d, 3) ls.add_to_constraints( X(3 * opp_h.to_corner() + d)-X(3 * h.from_corner() + d) );
		}
	}


	// assert symmetry
	for (auto f : m.iter_facets()) if (f.opposite().active()) for (auto h : f.iter_halfedges()) {
		auto opp_h = h.opposite_c();
		if (!opp_h.active()) continue;
		auto perm = permute_Jj_to_approx_Ji(ff[f.cell()], ff[opp_h.cell()]);
		if (!cut[f]) um_assert(perm.mid == 0);
		auto inv = permute_Jj_to_approx_Ji(ff[opp_h.cell()], ff[f.cell()]);
		um_assert(perm.inverse().mid == inv.mid);
	}

	if (paramCC.debug_seamless_across_cuts) {
		Trace::step("debug_continuity_across_cuts");
		for (auto f : m.iter_facets()) if (cut[f]) if (f.opposite().active()) for(auto h:f.iter_halfedges()) {
			auto opp_h = h.opposite_c();

			AxisPermutation perm;
			perm = permute_Jj_to_approx_Ji(ff[f.cell()], ff[opp_h.cell()]);
			FOR(d, 3) {
				LinExpr constraint = X(3 * h.to_corner() + d) - X(3 * h.from_corner() + d);
				FOR(dd, 3) {
					int match = int(perm.get_mat()[d][dd]);
					if (match == 0)  continue; // just to optimize
					constraint +=  match * X(3 * opp_h.to_corner() + dd) - match * X(3 * opp_h.from_corner() + dd);
				}
				ls.add_to_constraints(constraint);
			}
		}
	}
	// add equality on boundary

	if (paramCC.enforce_locked_cell) {
		Trace::step("align with locked cells");
		EdgeGraph eg(m);
		EdgeAttribute<int> fail(eg,-1);

		for (auto f : m.iter_facets()) if (locked_cells[f.cell()]) for (auto h : f.iter_halfedges()) FOR(d, 3) {
			LinExpr constraint = X(3 * h.to_corner() + d) - X(3 * h.from_corner() + d) + (U[h.from_corner()][d] - U[h.to_corner()][d]);
			ls.rb.reduce(constraint);
			if (constraint.size() < 2) {
				if (constraint.size() > 0) {
					plop(constraint);
					fail[eg.edge_from_halfedge(h)] = 0;
					//static int max_drop = 10;
					//if (max_drop -->0) drop(Triangle3(f), "fail");
				}
				continue;
			}
			ls.add_to_constraints(constraint);
		}
		Drop<PolyLine,EdgeAttribute<int> >(eg,fail).apply_wireframe("fail");
	}

	if (paramCC.align_with_boundary) {
		Trace::step("align with domain boundary");
		for (auto f : m.iter_facets()) {
			if (f.opposite().active() || locked_cells[f.cell()]) continue;
			Triangle3 tri = Triangle3(f);
			FOR(lv, 3) tri[lv] = ff[f.cell()] * tri[lv];
			vec3 n = tri.normal();
			if (constraint_type[f]==2) n = vec3(0,0,1);
			else if (constraint_type[f]==1 || constraint_type[f]==3)
				n = (std::abs(n[0]) >std::abs(n[1]) )?vec3(1,0,0):vec3(0,1,0);
			else continue;


			for( auto h: f.iter_halfedges())FOR(d, 3) {
				LinExpr constraint = n[d] * X(3 * h.from_corner() + d) - n[d] * X(3 * h.to_corner() + d);
				ls.rb.reduce(constraint);
				if (constraint.size() < 2) continue;
				ls.add_to_constraints(constraint);
			}
		}
	}
}


void CubeCover::add_int_constraints(CellFacetAttribute<bool>& cut, ConstrainedLeastSquares& ls, CellCornerAttribute<vec3>& U_brush, bool force_boundary) {


	Trace::Section sec("int constraints");

	// add integer coord for singularity graph uvs

	if (paramCC.GP_on_singularities) {
		Trace::step("force int on singu");

		for (auto h : m.iter_halfedges()) {
			if (!ff.edge_is_singular(Volume::Halfedge(m, h))) continue;
			int coord = ff.singular_edge_stable_coordinate(Volume::Halfedge(m, h));
			if (coord < 0 || coord>2) continue;
			FOR(d, 3) {
				if (d == coord) continue;
				LinExpr constraint = X(3 * h.from_corner() + d) -std::round(U[h.from_corner()][d]);
				ls.rb.reduce(constraint);

				if (constraint.size() == 0) continue;
				if (constraint.size() == 1 && constraint.front().index == ls.nfree) continue;
				ls.add_to_constraints(constraint);
			}
		}
	}


	if (paramCC.snapped_to_boundary) {
		// add integer coord for boundary
		if (force_boundary) {
			Trace::step("force int on boundary");
			for (auto f : m.iter_facets()) {
				if (f.opposite().active()) continue;
				int locked_branch = 0;
				vec3 n = Triangle3(f).normal();
				for (int d : { 1, 2 }) if (std::abs(ff[f.cell()].col(d) * n) > std::abs(ff[f.cell()].col(locked_branch) * n)) locked_branch = d;
				if (std::abs(ff[f.cell()][locked_branch] * n) < .9) continue;
				FOR(lh, 3) {
					auto h = f.halfedge(lh);
					FOR(d, 3) {
						LinExpr constraint = X(3 * h.from_corner() + locked_branch) + std::round(U[h.from_corner()][locked_branch]);
						ls.rb.reduce(constraint);
						if (constraint.size() < 2) continue;
						ls.add_to_constraints(constraint);
					}
				}
			}
		}
	}

	if (paramCC.debug_GP_across_cut) {
		for (auto h : m.iter_halfedges()) {
			auto opp_c = h.opposite_c();
			if (!opp_c.active()) continue;
			if (!cut[h.facet()]) continue;
			FOR(d, 3) {
				AxisPermutation perm;
				perm = permute_Jj_to_approx_Ji(ff[h.cell()], ff[opp_c.cell()]);

				LinExpr constraint = -1*X(3 * h.from_corner() + d);
				FOR(dd, 3) {
					int match = int(perm.get_mat()[d][dd]);
					if (match == 0) continue; // just to optimize
					constraint+=match * X(3 * opp_c.to_corner() + dd);
				}
				um_assert(constraint.size() == 2);
				double rhs = 0;
				FOR(term, 2) rhs += U[constraint[term].index / 3][constraint[term].index % 3] * constraint[term].value;
				constraint += -std::round(rhs) * X(ls.nfree);
				ls.rb.reduce(constraint);

				if (constraint.size() < 2) continue;
				//std::cerr << constraint << std::endl;
				ls.add_to_constraints(constraint);
			}
		}
	}
}

void CubeCover::integrate_field(ConstrainedLeastSquares& ls, EdgeGraph& eg, double edge_length) {
	Trace::step("Setup energy");
	for (auto h : m.iter_halfedges()) {
		int c_from = h.from_corner();
		int c_to = h.to_corner();

		vec3 rhs = -(1. / edge_length) *  (ff[h.cell()] * (h.to().pos() - h.from().pos()));
		FOR(d, 3) {
			ls.add_to_energy(UM::Linear::X(3 * c_to + d)-X(3 * c_from + d)+ rhs[d] * X(ls.nfree));
		}
	}
	Trace::step("Solve");
	ls.solve();
	Trace::step("Get result");
	FOR(c, m.ncorners()) FOR(d, 3)U[c][d] = ls.value(3 * c + d);
	
}






void CubeCover::compute_GP() {
	paramCC.debug_GP_across_cut = true;
	//paramCC.debug_seamless_across_cuts = false;

	Trace::Section sec("CubeCover");
	EdgeGraph eg(m);
	CellFacetAttribute<bool> cut(m, true);



	//CellCornerAttribute<vec3> U_backup(m);
	//for (auto c : m.iter_corners()) U_backup[c] = U[c];
	//auto show_compatibility = [&](std::string pre = "") {
	//	for (auto c : m.iter_cells()) {
	//		if (!locked_cells[c]) continue;
	//		AxisPermutation perm = permute_Jj_to_approx_Ji(uvw_to_jacobian(c, U), uvw_to_jacobian(c, U_backup));
	//		vec3 t = perm.get_mat() * U_backup[c.corner(0)] - U[c.corner(0)];
	//		FOR(lv, 4) {
	//			double diff = (perm.get_mat() * U_backup[c.corner(lv)] - U[c.corner(lv)] - t).norm();
	//			if (diff > 1e-12)
	//				std::cerr << pre << " :  diff = " << diff
	//				<< " diffvect = "
	//				<< perm.get_mat() * U_backup[c.corner(lv)] - U[c.corner(lv)] - t
	//				<< " U = " << U[c.corner(lv)] - U[c.corner(0)] << "  =>   rotated backup " << perm.get_mat() * (U_backup[c.corner(lv)] - U_backup[c.corner(0)]) << "\n";
	//		}
	//	}
	//	};


	//ff.show_cubes("init_cube");
	cut_graph_3D_and_brush_ff(cut);
	//ff.show_cubes("brushed_cube");

	CellCornerAttribute<vec3> U_brush(m);
	for (auto c : m.iter_corners()) U_brush[c] = U[c];
	Drop(m, U).apply("brushed U");











	//return;
	Drop(m, locked_cells)._skip_value(false).apply("locked");

	double edge_length = ToolBox(m).ave_edge_size();
	FOR(pass, 2) {
		plop(pass);
		Trace::step("Build reduction matrix");
		ConstrainedLeastSquares ls(3 * m.ncorners());
		add_real_constraints(cut, ls);
		if (pass == 1)
			add_int_constraints(cut, ls, U_brush, false);
		Trace::step("Integrate");
		integrate_field(ls, eg, edge_length);
		//show_compatibility(std::to_string(pass));

		
		if (!Trace::drop_mesh_is_active) continue;
		//Drop(m, U).apply_neg_det("negdet");
		//Drop(m, U).apply_iso("iso", 1);
		Drop(m, U).apply("U pass");
	}
}


void CubeCover::compute_seamless() {
	paramCC.debug_GP_across_cut = false;
	//paramCC.debug_seamless_across_cuts = false;

	Trace::Section sec("CubeCover_seamless");
	EdgeGraph eg(m);
	CellFacetAttribute<bool> cut(m, true);

	Drop(m, U).apply("inU");

	cut_graph_3D_and_brush_ff(cut);
	
	Drop(m, U).apply("brushU");

	//return;
	Drop(m, locked_cells)._skip_value(false).apply("locked");

	double edge_length = ToolBox(m).ave_edge_size();
	Trace::step("Build reduction matrix");
	ConstrainedLeastSquares ls(3 * m.ncorners());
	add_real_constraints(cut, ls);
	Trace::step("Integrate");
	integrate_field(ls, eg, edge_length);
	Drop(m, U).apply("seamlessU");
	if (!Trace::drop_mesh_is_active) return;

	Drop(m, U).apply_neg_det("negdet");
	Drop(m, U).apply_iso("iso", 1);
	Drop(m, U).apply("U pass");
}

