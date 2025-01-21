#include "seamless_to_GP.h"

#include <fullhex/gp_basic.h>

#include <framework/trace.h>

#include <toolbox.h>
#include <drop_attribute.h> 
#include <volume/frame3D.h>





using namespace UM::Linear;





SeamlessToGP::SeamlessToGP(Tetrahedra& p_m, CellCornerAttribute<vec3>& U, CellFacetAttribute<int>& constraint_type) : m(p_m), U(U), cut(m),constraint_type(constraint_type) {}

void SeamlessToGP::cut_graph_3D_and_brush_ff() {
	Trace::Section sec("cut_graph_3D");
	int seed = -1;
	Trace::step("compute dist 2 feature");
	CellAttribute<double> dist(m, 1e20);
	{
		auto cmp = [&](int left, int right) { return dist[left] > dist[right]; };
		std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);

		for (auto f : m.iter_facets()) if (!f.opposite().active() ) {
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
				GPTransitionFunction fct(f,U);
				
				FOR(lv, 4) U[opp_c.corner(lv)] = fct.apply(U[opp_c.corner(lv)]);

				cut[opp_f] = false;
				cut[f] = false;
				visited[opp_c] = true;
				queue.push(opp_c);
			}
		}
	}


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
			if (uvw_singu(h,U)) continue;
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
	}

}

void SeamlessToGP::add_real_constraints(ConstrainedLeastSquares& cls, bool force_boundary) {
	
	Trace::Section sec("real constraints");
	
		Trace::step("debug_continuity_outside_cuts");
		for (auto f : m.iter_facets()) if (!cut[f]) FOR(lh, 3) {
			auto h = f.halfedge(lh);
			auto opp_h = h.opposite_c();
			if (!opp_h.active()) continue;
			FOR(d, 3) cls.add_to_constraints( X(3 * opp_h.to_corner() + d)-X(3 * h.from_corner() + d) );
		}



		Trace::step("debug_continuity_across_cuts");
		for (auto f : m.iter_facets()) if (cut[f]) if (f.opposite().active()) for(auto h:f.iter_halfedges()) {
			auto opp_h = h.opposite_c();
			GPTransitionFunction fct(f,U);

			FOR(d, 3) {
				LinExpr constraint = X(3 * h.to_corner() + d) - X(3 * h.from_corner() + d);
				FOR(dd, 3) {
					int match = int(fct.ap.get_mat()[d][dd]);
					if (match == 0)  continue; // just to optimize
					constraint +=  match * X(3 * opp_h.to_corner() + dd) - match * X(3 * opp_h.from_corner() + dd);
				}
				cls.add_to_constraints(constraint);
			}
		}



		Trace::step("align with domain boundary");
		if (force_boundary) for (auto f : m.iter_facets()) {
			if (f.opposite().active()) continue;
			Triangle3 tri = uvw_tri(U,f);
			vec3 n = tri.normal();
			if (tri.cross_product().norm2() < 1e-10) continue;

			for( auto h: f.iter_halfedges())FOR(d, 3) {
				if (std::abs(n[d])<.9) continue;
				LinExpr constraint = n[d] * X(3 * h.from_corner() + d) - n[d] * X(3 * h.to_corner() + d);
				//ls.rb.leading_to_free(constraint);
				//if (constraint.size() < 2) continue;

				auto [already_satisfied,impossible] = constraint_status(cls,constraint);
				if (already_satisfied || impossible) continue;


				cls.add_to_constraints(constraint);
			}
		}
	
}


void SeamlessToGP::add_int_constraints(ConstrainedLeastSquares& cls,  bool force_boundary) {


	Trace::Section sec("int constraints");

	// add integer coord for singularity graph uvs


		Trace::step("force int on singu");

		for (auto h : m.iter_halfedges()) {
			if (!uvw_singu(Volume::Halfedge(m, h),U)) continue;

			FOR(d, 3) {
				if (std::abs(U[h.to_corner()][d]-U[h.from_corner()][d])>1e-5) continue;
				LinExpr constraint = X(3 * h.from_corner() + d) -std::round(U[h.from_corner()][d]);
				//cls.rb.leading_to_free(constraint);

				//if (constraint.size() == 0) continue;
				//if (constraint.size() == 1 && constraint.front().index == ls.nfree) continue;
				auto [already_satisfied,impossible] = constraint_status(cls,constraint);
				if (already_satisfied || impossible) continue;

				cls.add_to_constraints(constraint);
			}
		}
	


		// add integer coord for boundary
		if (force_boundary) {
			Trace::step("force int on boundary");
			for (auto f : m.iter_facets()) {
				if (f.opposite().active()) continue;
				vec3 n = uvw_tri(U,f).normal();
				
				FOR(branch,3) if (std::abs(n[branch])>.9) FOR(lh, 3) {
					auto h = f.halfedge(lh);
					FOR(d, 3) {
						LinExpr constraint = X(3 * h.from_corner() + branch) - std::round(U[h.from_corner()][branch]);
						cls.rb.leading_to_free(constraint);
						//if (constraint.size() < 1 || (constraint.size()==1 && constraint.front().index==-1)) continue;
						auto [already_satisfied,impossible] = constraint_status(cls,constraint);
						if (already_satisfied || impossible) continue;

						cls.add_to_constraints(constraint);
					}
				}
			}
		
	}

		for (auto h : m.iter_halfedges()) {
			auto opp_c = h.opposite_c();
			if (!opp_c.active()) continue;
			if (!cut[h.facet()]) continue;
			GPTransitionFunction fct(h.facet(),U);
			FOR(d, 3) {
				AxisPermutation perm = fct.ap;
				//perm = permute_Jj_to_approx_Ji(ff[h.cell()], ff[opp_c.cell()]);

				LinExpr constraint = -1*X(3 * h.from_corner() + d);
				FOR(dd, 3) {
					int match = int(perm.get_mat()[d][dd]);
					if (match == 0) continue; // just to optimize
					constraint+=match * X(3 * opp_c.to_corner() + dd);
				}
				um_assert(constraint.size() == 2);
				double rhs = 0;
				FOR(term, 2) rhs += U[constraint[term].index / 3][constraint[term].index % 3] * constraint[term].value;
				constraint = constraint-std::round(rhs);
				//cls.rb.leading_to_free(constraint);
				//if (constraint.size() < 2) continue;
				auto [already_satisfied,impossible] = constraint_status(cls,constraint);
				if (already_satisfied || impossible) continue;

				cls.add_to_constraints(constraint);
			}
		}
	
}

void SeamlessToGP::integrate_field(ConstrainedLeastSquares& ls, EdgeGraph& eg) {
	Trace::step("Setup energy");
	for (auto h : m.iter_halfedges()) {
		int c_from = h.from_corner();
		int c_to = h.to_corner();
		vec3 rhs = - (U[h.to_corner()]- U[h.from_corner()]);
		FOR(d, 3)  ls.add_to_energy(X(3 * c_to + d) - X(3*c_from+d) + rhs[d]);
		
	}
	Trace::step("Solve");
	ls.solve();
	Trace::step("Get result");
	FOR(c, m.ncorners()) FOR(d, 3)U[c][d] = ls.value(3 * c + d);
	
}






void SeamlessToGP::compute(bool force_boundary) {


			Drop(m, U).apply_iso("input iso 1");

	Trace::Section sec("SeamlessToGP");
	EdgeGraph eg(m);
	cut.fill(true);

	EdgeAttribute<bool> singu(eg);
	for(auto e:eg.iter_edges()) {auto h = eg.halfedge_from_edge(e);singu[e] = uvw_singu(h,U);}
	Drop<PolyLine,EdgeAttribute<bool> >(eg,singu)._skip_value(false).apply_wireframe("singu_uvw");

	cut_graph_3D_and_brush_ff();

	CellCornerAttribute<vec3> U_brush(m);
	for (auto c : m.iter_corners()) U_brush[c] = U[c];
	Drop(m, U).apply("brushed U");

	


	double edge_length = ToolBox(m).ave_edge_size();
	//FOR(pass, 2) {
		//plop(pass);
		Trace::step("Build reduction matrix");
		ConstrainedLeastSquares cls(3 * m.ncorners());
		add_real_constraints(cls, force_boundary);
		//if (pass == 1) 
		add_int_constraints(cls, false);


			for(auto c:m.iter_corners())  FOR(d,3) SparseVector v = express_with_free_variables(cls,X(3*c+d));

		Trace::step("Integrate");
		integrate_field(cls, eg);

		if (!Trace::drop_mesh_is_active) return;
		
		Drop(m, U).apply_neg_det();
		Drop(m, U).apply_iso("iso 1");

	//}
}
