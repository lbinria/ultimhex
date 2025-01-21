#include "FF_to_seamless.h"

#include <fullhex/gp_basic.h>

#include <framework/trace.h>

#include <toolbox.h>
#include <drop_attribute.h>
#include <volume/frame3D.h>





using namespace UM::Linear;





FFToSeamless::FFToSeamless(Tetrahedra& p_m, FF3D& p_ff, CellCornerAttribute<vec3>& U,CellFacetAttribute<int>& constraint_type) : m(p_m), ff(p_ff), U(U),constraint_type(constraint_type) {}

void FFToSeamless::cut_graph_3D_and_brush_ff(CellFacetAttribute<bool>& cut) {
	Trace::Section sec("cut_graph_3D");
	int seed = -1;
	Trace::step("compute dist 2 feature");
	CellAttribute<double> dist(m, 1e20);
	{
		auto cmp = [&](int left, int right) { return dist[left] > dist[right]; };
		std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);

		for (auto f : m.iter_facets()) if (constraint_type[f] ==5){
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

	
	for (auto c:m.iter_cells()) {
		mat3x3 Ju=uvw_to_jacobian(c, U);
		auto perm = permute_Jj_to_approx_Ji( ff[c],Ju);
		FOR(lv, 4) U[c.corner(lv)] = perm.get_mat() * U[c.corner(lv)];
	}

	ff.update_axis_permutation();

	Trace::step("remove boundaries and contractible loops ");

	for (auto f : m.iter_facets()) if (!f.opposite().active()) cut[f] = false;
	
	EdgeGraph eg(m);
	EdgeAttribute<int> nb_cuts(eg,0);
	EdgeAttribute<bool> on_border(eg,false);

	//precompute eg.edge_from_halfedge  (on gagne un facteur 2 sur la fonction)
	std::vector<int> h2e(m.nfacets()*3);
	for(auto h:m.iter_halfedges()) h2e[h]=eg.edge_from_halfedge(h);

	for(auto h:m.iter_halfedges()) if (cut[h.facet()])					nb_cuts [h2e[h]]++;
	for(auto h:m.iter_halfedges()) if (!h.facet().opposite().active())	on_border[h2e[h]]=true;
	bool done = false;
	while (!done) {
		done = true;
		for (auto f : m.iter_facets()) if (f.opposite().active() ) for(auto h:f.iter_halfedges()) if( cut[f]){
			if (nb_cuts[h2e[h]]!=1) continue; 
			if (on_border[h2e[h]]) continue; 
			if (ff.edge_is_singular(h)) continue;
			cut[f] = false;
			cut[f.opposite()] = false;
			for(auto cir : f.iter_halfedges())				nb_cuts[h2e[cir]]--;
			for(auto cir : f.opposite().iter_halfedges())	nb_cuts[h2e[cir]]--;
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
		Drop(m, U).apply("brushU");

	}

}





bool FFToSeamless::add_real_constraints(CellFacetAttribute<bool>& cut, ConstrainedLeastSquares& cls, bool force_boundary) {
	
	Trace::Section sec("real constraints");
	
		Trace::step("continuity_outside_cuts");
		for (auto f : m.iter_facets()) if (!cut[f] && f.opposite().active()) FOR(lh, 3) {
			auto h = f.halfedge(lh);
			auto opp_h = h.opposite_c();
			um_assert(opp_h.active()) ;
			FOR(d, 3) cls.add_to_constraints( X(3 * opp_h.to_corner() + d)-X(3 * h.from_corner() + d) );
		}

		//return true;
	// assert symmetry
	//for (auto f : m.iter_facets()) if (f.opposite().active()) for (auto h : f.iter_halfedges()) {
	//	auto opp_h = h.opposite_c();
	//	if (!opp_h.active()) continue;
	//	auto perm = permute_Jj_to_approx_Ji(ff[f.cell()], ff[opp_h.cell()]);
	//	if (!cut[f]) um_assert(perm.mid == 0);
	//	auto inv = permute_Jj_to_approx_Ji(ff[opp_h.cell()], ff[f.cell()]);
	//	um_assert(perm.inverse().mid == inv.mid);
	//}

		Trace::step("continuity_across_cuts");
		for (auto f : m.iter_facets()) if (cut[f]) if (f.opposite().active()) for(auto h:f.iter_halfedges()) {
			auto opp_h = h.opposite_c();

			AxisPermutation perm;
			perm = permute_Jj_to_approx_Ji(ff[f.cell()], ff[opp_h.cell()]);
			FOR(d, 3) {
				auto constraint = X(3 * h.to_corner() + d) - X(3 * h.from_corner() + d);
				FOR(dd, 3) {
					int match = int(perm.get_mat()[d][dd]);
					if (match == 0)  continue; // just to optimize
					constraint +=  match * X(3 * opp_h.to_corner() + dd) - match * X(3 * opp_h.from_corner() + dd);
				}
				cls.add_to_constraints(constraint);
			}
	}
	// add equality on boundary

		Trace::step("align with locked cells");
		EdgeGraph eg(m);
		EdgeAttribute<int> fail(eg,-1);int nfails = 0;

		for (auto f : m.iter_facets()) 
			if (constraint_type[f] ==-2)
				for (auto h : f.iter_halfedges()) FOR(d, 3) {
			LinExpr constraint = X(3 * h.to_corner() + d) - X(3 * h.from_corner() + d) + (U[h.from_corner()][d] - U[h.to_corner()][d]);
			auto [already_satisfied,impossible] = constraint_status(cls,constraint);
			if (impossible) {
					drop_triangle(Triangle3(f), "fail");
					//return false;
				
				continue;
			}
			cls.add_to_constraints(constraint);
		}
		if (nfails>0) Drop<PolyLine,EdgeAttribute<int> >(eg,fail).apply_wireframe("fail");

		Trace::step("align with domain boundary");
		if (force_boundary) for (auto f : m.iter_facets()) {
			if (f.opposite().active() ) continue;
			if (constraint_type[f]==-2) continue;
			Triangle3 tri = Triangle3(f);
			FOR(lv, 3) tri[lv] = ff[f.cell()] * tri[lv];
			vec3 n = tri.normal();
			if (constraint_type[f]==2) n = vec3(0,0,1);
			else if (constraint_type[f]==5)
				n = (std::abs(n[0]) >std::abs(n[1]) )?vec3(1,0,0):vec3(0,1,0);
			else continue;


			for( auto h: f.iter_halfedges())FOR(d, 3) {
				LinExpr constraint = n[d] * X(3 * h.from_corner() + d) - n[d] * X(3 * h.to_corner() + d);
				auto [already_satisfied,impossible] = constraint_status(cls,constraint);
				if (already_satisfied || impossible) continue;
				//ls.rb.leading_to_free(constraint);
				//if (constraint.size() < 2) continue;
				cls.add_to_constraints(constraint);
			}
		
	}
	return true;
}

bool FFToSeamless::integrate_field(ConstrainedLeastSquares& ls, double edge_length) {
	Trace::step("Setup energy");
	LinExpr expr;
	for (auto h : m.iter_halfedges()) {
		int c_from = h.from_corner();
		int c_to = h.to_corner();

		vec3 rhs = -(1. / edge_length) *  (ff[h.cell()] * (h.to().pos() - h.from().pos()));
		FOR(d, 3) ls.add_to_energy(X(3 * c_to + d)-X(3 * c_from + d)+rhs[d]);
	}
	Trace::step("Solve");
	ls.solve();
	Trace::step("Get result");
	FOR(c, m.ncorners()) FOR(d, 3)U[c][d] = ls.value(3 * c + d);
	return true;
}

bool FFToSeamless::compute(bool force_boundary) {

	Trace::Section sec("FFToSeamless_seamless");
	EdgeGraph eg(m);
	CellFacetAttribute<bool> cut(m, true);
	

	TRACE_ON(cut_graph_3D_and_brush_ff(cut));


	double edge_length = ToolBox(m).ave_edge_size();
	Trace::step("Build reduction matrix");
	ConstrainedLeastSquares cls(3 * m.ncorners());
	if (!add_real_constraints(cut, cls, force_boundary)) return false;


	//CellCornerAttribute<vec3> cid(m);
	//for(auto c:m.iter_corners()) {
	//	FOR(d,3){
	//	SparseVector v = express_with_free_variables(cls,X(3*c+d));
	//	um_assert(v.size()==1);
	//	cid[c][d] = v[0].index;
	//	um_assert(v[0].value==1);
	//	}
	//}
	//Drop(m,cid).apply("cid");
	//for(auto c:m.iter_corners())  FOR(d,3) SparseVector v = express_with_free_variables(cls,X(3*c+d));
	FOR(i,cls.rb.C.nrows()-1) {LinExpr le = X(i);cls.rb.leading_to_free(le);}

	Trace::step("Integrate");
	if (!integrate_field(cls, edge_length)) return false;

	if (Trace::drop_mesh_is_active){ 
		Drop(m, U).apply_neg_det("SEAMLESS_negdet");
		Drop(m, U).apply("SEAMLESS_U");
	}
	return true;
}

