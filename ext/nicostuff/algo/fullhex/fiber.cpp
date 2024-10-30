#include <fullhex/fiber.h>
#include <fullhex/gp_basic.h>
#include <framework/trace.h>

#include <toolbox.h>
#include <drop_attribute.h>


void permute_axis(AxisPermutation ap, std::array<bool, 3>& data) {
	std::array<bool, 3> nv={false,false,false};
	FOR(branch, 3) FOR(d, 3) if (ap.get_mat()[d][branch] != 0) nv[d] = data[branch];
	data = nv;
}

void crunch_selection_until_manifold(Tetrahedra& m, CellAttribute<bool>& selection) {

	EdgeGraph eg(m);
	bool done = false;
	while (!done) {
		done = true;
		Trace::step("Fix non manifold edges");
		for (auto e : eg.iter_edges()) {
			auto h = eg.halfedge_from_edge(e);
			int nb_bound = 0;
			for (auto cir : h.iter_CCW_around_edge()) {
				if (!selection[cir.cell()]) continue;
				auto opp_f = cir.facet().opposite();
				if (!opp_f.active()) nb_bound++;
				else if (!selection[opp_f.cell()]) nb_bound++;
			}
			if (nb_bound > 1) {
				done = false;
				for (auto cir : h.iter_CCW_around_edge()) if (!selection[cir.cell()]) std::cerr << ".";
				for (auto cir : h.iter_CCW_around_edge()) selection[cir.cell()] = false;
			}
		}
		if (!done) continue;



		Trace::step("Fix non manifold vertices");
		{
			DisjointSet ds(m.ncells());
			for (auto v : eg.iter_vertices()) {
				std::fill(ds.m_size.begin(), ds.m_size.end(), 1);
				std::iota(ds.m_ids.begin(), ds.m_ids.end(), 0);
				for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
					auto opp_f = cir.facet().opposite();
					if (!opp_f.active()) continue;
					if (selection[cir.cell()] == selection[opp_f.cell()]) ds.merge(opp_f.cell(), cir.cell());
				}
				int root_id[2] = {-1,-1};
					for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
						int side= selection[cir.cell()]? 0:1;
						int cell_root = ds.root(cir.cell());
						if (root_id[side] == -1) root_id[side] = cell_root;
						else if (root_id[side] != cell_root) root_id[side] = -2;
					}
					FOR(side,2)if (root_id[side] == -2) {
						done = false;
						for (auto e : v.iter_edges()) for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) selection[cir.cell()] = false;
					}
			}
		}
	}

	Tetrahedra sub_show;
	SubTetrahedra sub(m,sub_show, selection);
	TetBoundary bound(sub.proxy);
	DropSurface(bound.tri).apply("selection");
}



void keep_only_largest_connexe_component(Tetrahedra& m, CellAttribute<bool>& selection) {

	DisjointSet ds(m.ncells());
	for(auto f:m.iter_facets()) if (!f.on_boundary()){
		if (selection[f.cell()] && selection[f.opposite().cell()])
			ds.merge(f.cell(),f.opposite().cell());
	}
	std::vector<int> grp;
	int ngrp = ds.get_sets_id(grp);
	std::vector<int> grp_size(ngrp,0);
	for(auto c:m.iter_cells()) grp_size[grp[c]]++;
	int largest_grp=0;
	FOR(i,ngrp) if (grp_size[largest_grp]<grp_size[i]) largest_grp = i;

	for(auto c:m.iter_cells()) selection[c] = (grp[c] ==largest_grp);
	Tetrahedra subshow;
	SubTetrahedra sub(m,subshow, selection);
	TetBoundary bound(sub.proxy);
	DropSurface(bound.tri).apply("largestCC");
}

void grow_selection_until_manifold(Tetrahedra& m, CellAttribute<bool>& selection) {

	{
		Trace::step("Dilate selection");
		CellAttribute<bool> tmp(m);
		FOR(it, 5) {
			for (auto c : m.iter_cells()) tmp[c] = selection[c];
			for (auto f : m.iter_facets())
				if (selection[f.cell()] && f.opposite().active())
					tmp[f.opposite().cell()] = true;
			for (auto c : m.iter_cells()) selection[c] = tmp[c];
		}
	}
	EdgeGraph eg(m);
	bool done = false;
	while (!done) {
		done = true;

		Trace::step("Fix non manifold edges");
		for (auto e : eg.iter_edges()) {
			auto h = eg.halfedge_from_edge(e);
			int nb_bound = 0;
			for (auto cir : h.iter_CCW_around_edge()) {
				if (!selection[cir.cell()]) continue;
				auto opp_f = cir.facet().opposite();
				if (!opp_f.active()) nb_bound++;
				else if (!selection[opp_f.cell()]) nb_bound++;
			}
			if (nb_bound > 1) {
				done = false;
				for (auto cir : h.iter_CCW_around_edge()) if (!selection[cir.cell()]) std::cerr << ".";
				for (auto cir : h.iter_CCW_around_edge()) selection[cir.cell()] = true;
			}
		}
		if (!done) continue;



		Trace::step("Fix non manifold vertices");
		{
			DisjointSet ds(m.ncells());
			DisjointSet ds_edge(eg.nedges());
			for (auto v : eg.iter_vertices()) {
				std::fill(ds.m_size.begin(), ds.m_size.end(), 1);
				std::iota(ds.m_ids.begin(), ds.m_ids.end(), 0);
				for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
					auto opp_f = cir.facet().opposite();
					if (!opp_f.active()) continue;
					if (!selection[cir.cell()]) continue;
					if (selection[opp_f.cell()]) ds.merge(opp_f.cell(), cir.cell());
				}


				{// check only one connex component
					int root_id = -1;
					for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
						if (!selection[cir.cell()]) continue;
						int cell_root = ds.root(cir.cell());
						if (root_id == -1) root_id = cell_root;
						else if (root_id != cell_root) root_id = -2;
					}
					if (root_id == -2) {
						done = false;
						for (auto e : v.iter_edges()) for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) selection[cir.cell()] = true;
						continue;
					}
					if (root_id == -1) {
						continue;
					}
				}

				std::fill(ds_edge.m_size.begin(), ds_edge.m_size.end(), 1);
				std::iota(ds_edge.m_ids.begin(), ds_edge.m_ids.end(), 0);

				for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
					auto opp_f = cir.facet().opposite();
					if (!selection[cir.cell()]) continue;
					if (opp_f.active() && selection[opp_f.cell()]) continue;
					ds_edge.merge(e, eg.edge_from_halfedge(cir.prev().opposite_f()));
				}
				std::set<int> edge_roots;
				for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
					auto opp_f = cir.facet().opposite();
					if (!selection[cir.cell()]) continue;
					if (opp_f.active() && selection[opp_f.cell()]) continue;
					edge_roots.insert(ds_edge.root(e));
				}

				if (edge_roots.size() > 1) {
					plop(edge_roots.size());
					done = false;
					for (auto e : v.iter_edges()) for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) selection[cir.cell()] = true;
					continue;
				}
			}
		}
	}
	Tetrahedra subshow;
	SubTetrahedra sub(m, subshow,selection);
	TetBoundary bound(sub.proxy);
	DropSurface(bound.tri).apply("selection");
}

struct MergeSets {
	MergeSets(int nb_elts) { next.resize(nb_elts); FOR(i, nb_elts) next[i] = i; }
	bool same_set(int i, int j) {
		int it = i;
		do {
			if (it == j) return true;
			it = next[it];
		} while (it != i);
		return false;
	}
	void merge(int i, int j) {
		if (!same_set(i, j)) std::swap(next[i], next[j]);
	}
	std::vector<int> next;
};

void show_map_singu(Tetrahedra& m, CellCornerAttribute<vec3>& U) {
	EdgeGraph eg(m);
	{// check uvw_sing works
		EdgeAttribute<bool> singu(eg, 0);
		for (auto e : eg.iter_edges()) {
			auto h = eg.halfedge_from_edge(e);
			singu[e] = uvw_singu(h, U);
		}
		Drop<PolyLine, EdgeAttribute<bool> >(eg, singu)._skip_value(false).apply_wireframe("uvw_sing");
	}
	CellFacetAttribute<int> ap(m, 0);
	for (auto f : m.iter_facets())if (f.opposite().active()) {
		GPTransitionFunction tf(f, U);
		ap[f] = tf.ap.mid;
	}
	Drop(m, ap)._skip_value(0).apply("ap");
}

void FiberSegmentation::trace_fiber() {

	CellAttribute<bool> touched(m);
	int num_fiber = 0;

	for (auto c : m.iter_cells()) fiber[c] = selection[c] ? -1 : -2;


	CellFacetAttribute<int> fiber_cross(m,-1);

	for (auto f_seed : m.iter_facets()) {
		if (f_seed.opposite().active()) continue;
		auto c = f_seed.cell();
		if (!selection[c]) continue;
		//plop(num_fiber);

		if (constraint_type[f_seed]!=2 || uvw_tri(U, f_seed).normal()[2] > 0) continue;
		fiber_cross[f_seed] = 0;
		touched.fill(false);


		BBox2 bbox_f_seed;
		for (auto h : f_seed.iter_halfedges()) bbox_f_seed.add(U[h.from_corner()].xy());
		Triangle2 tri_f_seed = uvw_tri(U, f_seed).xy();

		std::vector<int> new_fiber_elt;
		auto add_to_fiber_if_in_prism = [&](Volume::Cell cell_to_test) {
			if (fiber[cell_to_test] != -1) return false;
			vec3 bc = tri_f_seed.bary_coords(uvw_tet(U, cell_to_test).bary_verts().xy());
			if (bc[0] > 0 && bc[1] > 0 && bc[2] > 0) {
				new_fiber_elt.push_back(cell_to_test);
				return true;
			} 
			return false;
			};


		bool num_fiber_was_used = false;
		std::vector<int> cells_to_visit;
		cells_to_visit.push_back(c);
		touched[c] = true;
		num_fiber_was_used = add_to_fiber_if_in_prism(c) || num_fiber_was_used;

		bool traverse_mesh = true;
		double max_cross = -1e20;
		while (!cells_to_visit.empty()) {
			c = cells_to_visit.back();
			cells_to_visit.pop_back();
			for (auto f : c.iter_facets()) {
				auto opp_f = f.opposite();
				if (!opp_f.active()) {
					//plop(constraint_type[f]);plop(uvw_tet(U,f.cell()).volume());plop(uvw_tri(U, f).normal()[2]);
					bool crossed =  (constraint_type[f]==2 && uvw_tet(U,f.cell()).volume()>0 && uvw_tri(U, f).normal()[2] > 0);
				//	plop(crossed);
					if (crossed) {
						fiber_cross[f] = 1;
						max_cross = std::max(max_cross,uvw_tri(U,f).bary_verts()[2]);
					}
					//traverse_mesh = traverse_mesh || crossed;					
					continue;
				}
				auto opp_c = opp_f.cell();
				if (touched[opp_c]) continue;


				GPTransitionFunction tf(f, U);
				if (!tf.ap.is_identity()) { traverse_mesh = false; break; }
				FOR(lv, 4) U[opp_c.corner(lv)] = tf.apply(U[opp_c.corner(lv)]);

				BBox2 bbox_opp_c;
				FOR(lv, 4)  bbox_opp_c.add(U[opp_c.corner(lv)].xy());
				if (bbox_f_seed.intersect(bbox_opp_c)) {
					touched[opp_c] = true;
					cells_to_visit.push_back(opp_c);
					num_fiber_was_used = add_to_fiber_if_in_prism(opp_c) || num_fiber_was_used;
				}
			}
		}
		for (int nc : new_fiber_elt) 
			traverse_mesh = traverse_mesh && max_cross > uvw_tet(U,Volume::Cell(m,nc)).bary_verts()[2];
		
		if (num_fiber_was_used && traverse_mesh) {
			for (int nc : new_fiber_elt) fiber[nc] = num_fiber;
			num_fiber++;
		}
	}

	for (auto c : m.iter_cells()) if (fiber[c] ==-1) {
		for (auto f:c.iter_facets()) if (!f.on_boundary() && fiber[f.opposite().cell()]>=0)fiber[c] =-3;
	}

	for (auto c : m.iter_cells()) {
		if (fiber[c] == -2)fiber[c] = -1;
		if (fiber[c] == -3) fiber[c] = num_fiber++;
	}
	Drop(m, fiber)._skip_value(-1).apply("traced fibers");
	Drop(m,fiber_cross).apply("fibercross");
}

// the new version brush the U
void FiberSegmentation::compute_singu_constraints(CellAttribute<std::array<bool, 3> >& constraint) {

	for (auto c : m.iter_cells()) FOR(d, 3) constraint[c][d] = false;

	EdgeGraph eg(m);
	EdgeAttribute<bool> edge_singu(eg);
	for (auto e : eg.iter_edges()) edge_singu[e] = uvw_singu(eg.halfedge_from_edge(e), U);

	for (auto e : eg.iter_edges()) if (edge_singu[e]) for (auto seed_h : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
		auto seed_c = seed_h.cell();

		{
			// rotate to put stable direction in Z, if avalaible
			int stable_branch = -1;
			vec3 geom = U[seed_h.to_corner()] - U[seed_h.from_corner()];
			FOR(d, 3) if (std::abs(geom[(d + 1) % 3]) < 1e-20 && std::abs(geom[(d + 2) % 3]) < 1e-20) stable_branch = d;
			if (stable_branch == -1) continue;
			AxisPermutation ap;
			if (stable_branch == 0) ap.mid = 9;
			if (stable_branch == 1) ap.mid = 5;
			FOR(lv, 4) U[seed_c.corner(lv)] = ap.get_mat() * U[seed_c.corner(lv)];
			permute_axis(ap, constraint[seed_c]);
		}

		constraint[seed_c][2] = true;


		FOR(branch, 2) {
			vec3 P_uvw = .5 * (U[seed_h.to_corner()] + U[seed_h.from_corner()]);
			int branch_n = (branch + 1) % 2;


			if (constraint[seed_c][branch]) continue;



			CellAttribute<std::array<bool, 3> > local_constraint(m);
			for (auto c : m.iter_cells()) FOR(d, 3) local_constraint[c][d] = false;

			BBox1 range;
			FOR(lv, 4) range.add(vec<1>{U[seed_h.cell().corner(lv)][2]});


			std::vector<int> stack;
			stack.push_back(seed_c);

			while (!stack.empty()) {
				int cur = stack.back();
				stack.pop_back();

				if (local_constraint[cur][branch]) continue;
				local_constraint[cur][branch] = true;


				//if (DEBUG_max_iter--<0) break;// to see only init conditions
				vec3 n_uvw = mat3x3::identity()[branch_n];
				Volume::Cell c(m, cur);
				for (auto f : c.iter_facets()) {
					auto opp = f.opposite();
					if (!opp.active()) continue;


					// lets brush
					GPTransitionFunction fct(f, U);
					FOR(lv, 4) U[opp.cell().corner(lv)] = fct.apply(U[opp.cell().corner(lv)]);
					permute_axis(fct.ap, constraint[opp.cell()]);

					// endof lets brush



					vec3 bc(0, 0, 0);
					Triangle3 tr(U[f.halfedge(0).from_corner()], U[f.halfedge(1).from_corner()], U[f.halfedge(2).from_corner()]);
					Triangle3 tr_opp(U[f.halfedge(0).opposite_c().to_corner()], U[f.halfedge(1).opposite_c().to_corner()], U[f.halfedge(2).opposite_c().to_corner()]);

					vec3 val; FOR(d, 3) val[d] = tr[d][branch_n] - P_uvw[branch_n];

					if (val[0] * val[1] > 0 && val[1] * val[2] > 0 && val[2] * val[1] > 0) continue;

					stack.push_back(opp.cell());

				}
			}

			for (auto c : m.iter_cells()) {
				if (!local_constraint[c][0] && !local_constraint[c][1]) continue;
				for (auto f : c.iter_facets()) for (auto h : f.iter_halfedges()) {
					//auto e = eg.edge_from_halfedge(h);
					BBox1 box;
					box.add(vec<1>{U[h.from_corner()][2]});
					box.add(vec<1>{U[h.to_corner()][2]});
					if (edge_singu[e]) range.add(box);
				}
			}

			for (auto c : m.iter_cells()) {
				if (!local_constraint[c][0] && !local_constraint[c][1]) continue;
				BBox1 box;
				FOR(lv, 4) box.add(vec<1>{U[c.corner(lv)][2]});
				if (range.intersect(box))
					constraint[c][branch] = true;
			}
			//{
			//	double ave = ToolBox(m).ave_edge_size(); FOR(branch, 3) {
			//		CellAttribute<vec3> tmp(m, vec3(0, 0, 0)); for (auto c : m.iter_cells())if (constraint[c][branch]) tmp[c] = uvw_to_jacobian(c, U)[branch];
			//		Drop(m, tmp)._force_length(.5 * ave, false)._force_radius(.05 * ave, false).apply_arrow("constraint");
			//	}return;
			//}
		}
	}


	Trace::step("render");

	double ave = ToolBox(m).ave_edge_size();
	FOR(branch, 3) {
		CellAttribute<vec3> tmp(m, vec3(0, 0, 0));
		for (auto c : m.iter_cells())if (constraint[c][branch]) tmp[c] = uvw_to_jacobian(c, U)[branch];
		Drop(m, tmp)._force_length(.5 * ave, false)._force_radius(.05 * ave, false).apply_arrow("constraint");
	}
	//Drop(m,U).apply_iso("iso_output");

}

bool z_can_be_normal(vec3 n) {
	return std::abs(n[2]) > .9;
}
bool z_can_be_stable(vec3 n) {
	FOR(d, 2) if (std::abs(n[d]) < .1) return false;
	return std::abs(n[2]) < .2;
}
bool FiberSegmentation::orient_seed(
	CellAttribute<std::array<bool, 3> >& constraint, CellAttribute<std::array<bool, 3> >& in_previous,
	Volume::Facet seed_f_tet
) {
	auto c = seed_f_tet.cell();
	vec3 n = uvw_tri(U, seed_f_tet).normal();
	FOR(test_stable, 3) {
		mat3x3 rot = mat3x3::identity();
		if (test_stable == 0) rot = roty(M_PI / 2.);
		if (test_stable == 1) rot = rotx(M_PI / 2.);
		if (!constraint[c][test_stable] && !in_previous[c][test_stable] && z_can_be_stable(rot * n)) {
			if (test_stable != 2) {
				AxisPermutation ap;
				if (test_stable == 0) ap.mid = 9;
				if (test_stable == 1) ap.mid = 5;
				plop(test_stable);
				FOR(lv, 4) U[c.corner(lv)] = ap.get_mat() * (U[c.corner(lv)]);
				permute_axis(ap, constraint[c]);
				permute_axis(ap, in_previous[c]);
			}
			return true;
		}
	}
	return false;
}

bool FiberSegmentation::grow_from_seed(
	CellAttribute<std::array<bool, 3> >& constraint,
	CellAttribute<std::array<bool, 3> >& in_previous,
	Volume::Facet seed_f_tet,
	CellAttribute<std::array<bool, 3> >& in_current,
	CellFacetAttribute<int>& local_constraint_type,
	FacetAttribute<int>& chart_id, FacetAttribute<vec3>& stable_direction, int& nb_charts) {


	local_constraint_type.fill(-1);

	if (!orient_seed(constraint, in_previous, seed_f_tet)) return false;

	
	um_assert(!in_previous[seed_f_tet.cell()][2]);

	for (auto c : m.iter_cells()) FOR(d, 3) in_current[c][d] = false;

	CellAttribute<double> dist(m, 1e20);
	CellAttribute<vec2> seed_pos(m);

	auto cmp = [&](int left, int right) { return dist[left] > dist[right]; };
	std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);


	//int nb_fixables = 0;
	dist[seed_f_tet.cell()] = 0;
	seed_pos[seed_f_tet.cell()] = uvw_tet(U, seed_f_tet.cell()).bary_verts().xy();
	queue.push(seed_f_tet.cell());


	//grow the region
	Triangles seeds;
	while (!queue.empty()) {
		//std::cerr<<".";
		Volume::Cell c(m, queue.top());
		queue.pop();
		

		if (in_current[c][0]) continue;
		if (in_current[c][1]) continue;
		if (constraint[c][2]) continue;
		in_current[c][2] = true;

		for (auto f : c.iter_facets()) {
			auto opp_f = f.opposite();
			if (!opp_f.active()) {
				if (dist[c] > 0 && z_can_be_stable(uvw_tri(U, f).normal())) {
					ToolBox(seeds).add_triangle(Triangle3(f));
					dist[c] = 0;
					seed_pos[c] = uvw_tet(U, c).bary_verts().xy();
					queue.push(c);
					//plop(c);
				}
				continue;
			}
			GPTransitionFunction tf(f, U);
			if (dist[opp_f.cell()] < 1e20 && !tf.ap.is_identity()) continue;

			auto opp_c = opp_f.cell();
			FOR(lv, 4) U[opp_c.corner(lv)] = tf.apply(U[opp_c.corner(lv)]);

			permute_axis(tf.ap, constraint[opp_c]);
			permute_axis(tf.ap, in_current[opp_c]);
			permute_axis(tf.ap, in_previous[opp_c]);


			vec2 c_U = uvw_tet(U, c).bary_verts().xy();
			vec2 opp_c_U = uvw_tet(U, opp_c).bary_verts().xy();
			double nv_dist_acc = dist[c] + (c_U - opp_c_U).norm();
			double nv_dist_direct = (seed_pos[c] - opp_c_U).norm();
			double nv_dist =(nv_dist_acc>5) ? nv_dist_acc : nv_dist_direct;
			if (dist[opp_c] > nv_dist + 1e-10) {
				seed_pos[opp_c] = seed_pos[c];
				dist[opp_c] = nv_dist;
				//static int printstart=0; if (printstart++>1000000){plop(opp_c);plop(dist[opp_c]);}
				queue.push(opp_c);
			}
		}
	}
	DropSurface(seeds).apply("seeds");
	Drop(m, U).apply("U");
	for (auto c : m.iter_cells()) FOR(branch, 3) in_previous[c][branch] = in_previous[c][branch] || in_current[c][branch];

	for (auto c : m.iter_cells()) FOR(branch, 3) selection[c] = selection[c] || in_current[c][branch];


	Drop(m, selection)._skip_value(false).apply("selection");

 
	for (auto f : m.iter_facets()) {
		if (!in_current[f.cell()][2]) continue;
		if (f.opposite().active()) continue;
		if (std::abs(uvw_tri(U, f).normal()[2]) < .1) local_constraint_type[f] = 5;
		if (z_can_be_normal(uvw_tri(U, f).normal())) local_constraint_type[f] = 2;
	}
 	std::vector<bool> compatible_charts(nb_charts, true);
	for (auto f : bound.tri.iter_facets()) {
		auto f_tet = bound.tet_facet(f);
		if (local_constraint_type[f_tet] == 2 && chart_id[f] != -1)
			compatible_charts[chart_id[f]] = false;
	}
 	for (auto f : bound.tri.iter_facets()) {
		auto f_tet = bound.tet_facet(f);
		if (chart_id[f] != -1 && !compatible_charts[chart_id[f]]) {
			local_constraint_type[f_tet] = -2;
		}
	}


	return true;
}

void FiberSegmentation::compute_chart_id(FacetAttribute<int>& chart_id, FacetAttribute<vec3>& stable_direction, int& nb_charts) {

	FacetAttribute<bool> can_be_constraint(bound.tri, false);

	for (auto f : bound.tri.iter_facets()) {
		auto f_tet = bound.tet_facet(f);
		stable_direction[f] = vec3(0, 0, 0);
		vec3 n = uvw_tri(U, f_tet).normal();
		mat3x3 J = uvw_to_jacobian(f_tet.cell(), U);
		FOR(axe, 3) {

			vec3 rn = n;
			if (axe == 0) rn = roty(M_PI / 2.) * n;
			if (axe == 1) rn = rotx(M_PI / 2.) * n;

			if (z_can_be_stable(rn))		can_be_constraint[f] = true;
			else continue;

			stable_direction[f] = (J.invert() * mat3x3::identity()[axe]).normalized();
		}
	}
	Drop(bound.tri, stable_direction).apply_arrow("chart_dir");
	Drop(bound.tri, can_be_constraint).apply("can_be_constraint");

	DisjointSetWithNeutral dsn(bound.tri.nfacets());
	for (auto f : bound.tri.iter_facets()) if (stable_direction[f].norm2() == 0) dsn.set_neutral(f);

	for (auto h : bound.tri.iter_halfedges()) {
		auto f = h.facet();
		auto o = h.opposite().facet();
		if (!can_be_constraint[f] ) continue;
		if (!can_be_constraint[o] ) continue;
		if (std::abs(stable_direction[f] * stable_direction[o]) > .9 ) dsn.merge(f, o);
	}
	nb_charts = dsn.get_sets_id(chart_id.ptr->data);

}

void FiberSegmentation::grow_regions(CellAttribute<std::array<bool, 3> >& constraint) {


	FacetAttribute<int> chart_id(bound.tri);
	FacetAttribute<vec3> stable_direction(bound.tri, vec3(0, 0, 0));
	int nb_charts;
	compute_chart_id(chart_id, stable_direction, nb_charts);

	//DropVolume(m).apply("tet");
	//DropSurface(bound.tri).apply("surface");
	Drop(m,U).apply_iso("input U");
	Drop(bound.tri, chart_id).apply("charts");




	double best_score = 0;
	Volume::Facet best_seed(m, -1);

	CellAttribute<std::array<bool, 3> > in_previous(m);
	copy_attribute(constraint, in_previous);
	CellAttribute<std::array<bool, 3> > in_current(m);


		
	CellFacetAttribute<int> local_constraint_type(m);
	for (auto seed_f_tri : bound.tri.iter_facets()) {
		auto seed_f_tet = bound.tet_facet(seed_f_tri);

		//init the region
		if (!grow_from_seed(constraint, in_previous, seed_f_tet, in_current, local_constraint_type, chart_id, stable_direction, nb_charts)) continue;

		double score = 0;
		for (auto f : m.iter_facets()) {
			if (!in_current[f.cell()][2]) continue;
			if (f.opposite().active()) continue;
			if (z_can_be_stable(uvw_tri(U, f).normal())) score += 1;
		}
		plop(score);
		if (best_score < score) {
			best_score = score;
			best_seed = seed_f_tet;
		}
	}

	copy_attribute(constraint, in_previous);
	if (!best_seed.active()) return;

	grow_from_seed(constraint, in_previous, best_seed, in_current, constraint_type, chart_id, stable_direction, nb_charts);
	for (auto c : m.iter_cells()) selection[c] = in_current[c][2];

	Drop(m, U).apply("Selected U");
	Drop(m, selection)._skip_value(false).apply("Selected  selection");
	Drop(m, constraint_type)._skip_value(-2).apply("Selected  constraints type");

}


void FiberSegmentation::cleanup_fibers(){
	for (auto c : m.iter_cells()) selection[c] = fiber[c]>=0 ;
	Drop(m,selection)._skip_value(false).apply("fiber selection");
	crunch_selection_until_manifold(m, selection);
	Drop(m,selection)._skip_value(false).apply("crunched selection");
	keep_only_largest_connexe_component(m, selection);
	Drop(m,selection)._skip_value(false).apply("largest connex selection");
	
	for (auto c : m.iter_cells()) if (!selection[c])  fiber[c]=-1; 


	for (auto f : m.iter_facets()){
		auto opp = f.opposite();
		if (!opp.active() && fiber[f.cell()] == -1) constraint_type[f] = -1;
		if (opp.active() && fiber[f.cell()]>=0 && fiber[opp.cell()]<0)    constraint_type[f] = -2;
	}
}

bool FiberSegmentation::apply() {


	CellAttribute<std::array<bool, 3> > constraint(m, { false,false,false });
	m.connect();
	TRACE_OFF(compute_singu_constraints(constraint));

	TRACE_OFF(grow_regions(constraint));
	
	TRACE_OFF(trace_fiber());


	TRACE_ON(cleanup_fibers());

	//Drop(m,selection)._skip_value(false).apply("Final selection");
	Drop(m,fiber)._skip_value(-1).apply("SELECT_Fibers");
	Drop(m,constraint_type).apply("SELECT_Constraint type");
	

	// check if the selection is large enough
	// and if there is enough facets that are likely to be better aligned
	int select_size=0;for(auto c:m.iter_cells()) if (selection[c]) select_size++;
	int fixable_size=0;for(auto f:m.iter_facets()) if (constraint_type[f]==5) fixable_size++;
	return select_size>200 && fixable_size>50;


}bool FiberSegmentation::apply_to_deformed_only() {


	
	TRACE_ON(trace_fiber());


	//TRACE_ON(cleanup_fibers());

	Drop(m,fiber)._skip_value(-1).apply("SELECT_Fibers");
	Drop(m,constraint_type).apply("SELECT_Constraint type");
	return true;
}

