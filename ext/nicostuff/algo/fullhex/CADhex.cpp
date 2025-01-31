#include <fullhex/CADhex.h>
#include <fullhex/gp_basic.h>
#include <volume/hex_select.h>
#include<algo/dynamicKNN.h>

HexCAD::HexCAD(Hexahedra& hex) :hex(hex), emb(hex.points) { verbose_level(1); }


HexCAD& HexCAD::drop_SJ(std::string name) {
	if (verbose > 0) ToolBox(hex).drop_scaled_jacobien(name);
	return *this;
}


const vec3 hex_face_to_normal[6] = {
	vec3(-1, 0, 0), vec3(1, 0, 0),
	vec3(0, -1, 0),vec3(0, 1, 0),
	vec3(0, 0, -1),vec3(0, 0, 1)
};


using namespace Linear;


void generate_cubes(
	Tetrahedra& tet,
	CellCornerAttribute<vec3>& uvw,
	Hexahedra& hex,
	CellAttribute<int>& cube_tet_id,
	CellAttribute<vec3>& cube_pos_uvw,
	CellAttribute<vec3>& cube_pos,
	DynamicKNN<3>& cube_knn,
	double hex_size
) {

	for (auto tetra : tet.iter_cells()) {
		BBox3 box;
		Tetrahedron tet_uvw(uvw[tetra.corner(0)], uvw[tetra.corner(1)], uvw[tetra.corner(2)], uvw[tetra.corner(3)]);

		if (tet_uvw.volume() < 1e-10) continue;
		FOR(lv, 4) box.add(tet_uvw[lv]);
		box.dilate(1e-20);
		for (int i = int(std::floor(box.min[0])); i < int(std::ceil(box.max[0])); i++)
			for (int j = int(std::floor(box.min[1])); j < int(std::ceil(box.max[1])); j++)
				for (int k = int(std::floor(box.min[2])); k < int(std::ceil(box.max[2])); k++) {
					// test is_in in texture space
					vec3 guess_uvw = vec3(i, j, k) + .5 * vec3(1, 1, 1);
					vec4 bc = tet_uvw.bary_coords(guess_uvw);
					bool is_in = true;
					FOR(lv, 4) is_in = is_in && (bc[lv] > -1e-20);
					if (!is_in) continue;

					// create the hex centered on guess_uvw
					vec3 bary_pos(0, 0, 0);
					FOR(lv, 4) bary_pos += bc[lv] * tetra.vertex(lv).pos();
					bool already_in = false;


					auto nn = cube_knn.query(bary_pos);
					if (nn.size() > 0) {
						int check_c = nn[0];
						//FOR(check_c, hex.ncells()) 
						already_in = already_in || (cube_pos[check_c] - bary_pos).norm() < hex_size / 100.;
					}
					if (already_in) continue;

					int offset_v = hex.points.create_points(8);
					int offset_c = hex.create_cells(1);
					cube_pos_uvw[offset_c] = guess_uvw;
					cube_pos[offset_c] = bary_pos;
					cube_tet_id[offset_c] = tetra;
					FOR(lv, 8) hex.vert(offset_c, lv) = offset_v + lv;

					mat3x3 J_inv = uvw_to_jacobian(tetra, uvw).invert();
					FOR(di, 2)FOR(dj, 2)FOR(dk, 2)
						hex.points[offset_v + di + 2 * dj + 4 * dk]
						= cube_pos[offset_c] + .5 * J_inv * (vec3(di, dj, dk) - .5 * vec3(1, 1, 1));
				}
	}
}


void link_cubes(
	Tetrahedra& tet,
	CellCornerAttribute<vec3>& uvw,
	Hexahedra& hex,
	CellAttribute<int>& cube_tet_id,
	CellAttribute<vec3>& cube_pos_uvw,
	CellAttribute<vec3>& cube_pos,
	DynamicKNN<3>& cube_knn,
	CellFacetAttribute<int>& hex_facet_emb,
	double hex_size
) {
	DisjointSet ds(hex.nverts());
	hex.connect();
	Trace::step("Find links between cubes");
	// connect cube faces to "opposite OR surface point"
	FOR(cur_hex, hex.ncells()) FOR(lf, 6) {
		vec3 dir = hex_face_to_normal[lf];
		vec3 objective_uvw = cube_pos_uvw[cur_hex] + dir;
		Volume::Cell cur_tet(tet, cube_tet_id[cur_hex]);

		int objective_hex = -1;
		vec3 I;
		vec3 G = cube_pos_uvw[cur_hex];
		GPTransitionFunction tf;
		while (true) {
			Tetrahedron tet_uvw(
				tf.apply(uvw[cur_tet.corner(0)]),
				tf.apply(uvw[cur_tet.corner(1)]),
				tf.apply(uvw[cur_tet.corner(2)]),
				tf.apply(uvw[cur_tet.corner(3)])
			);
			if (tet_uvw.volume() < 0) {
				Trace::alert("Cannot cross a tet with negative volume");
				break;
			}

			vec4 bc = tet_uvw.bary_coords(objective_uvw);

			{// check if the destination is reached
				bool is_in = true;
				FOR(lv, 4) is_in = is_in && (bc[lv] > -1e-10);
				if (is_in) {
					vec3 objective_xyz = vec3(0, 0, 0);
					FOR(lv, 4) objective_xyz += bc[lv] * cur_tet.vertex(lv).pos();
					int inner_cur_hex = cube_knn.query(objective_xyz)[0];
					if ((cube_pos[inner_cur_hex] - objective_xyz).norm() < 1e-2 * hex_size && cube_tet_id[inner_cur_hex] == cur_tet)
						objective_hex = inner_cur_hex;
					else std::cerr << "Found a tet that contains the objective, but does not match an hex :( \n\n";
					break;
				}
			}


			{// objective_uvw is not in cut_tet, need to visit the next one.
				Volume::Facet out_facet(tet, -1);

				for (auto f : cur_tet.iter_facets()) {
					Triangle3 tr_uvw({
						tf.apply(uvw[f.corner(0)]),
						tf.apply(uvw[f.corner(1)]),
						tf.apply(uvw[f.corner(2)])
						});
					vec3 bc_tri;

					if (!Intersect::triangle_line(tr_uvw, G, dir, bc_tri)) continue;
					if (Tetrahedron(tr_uvw[0], tr_uvw[1], tr_uvw[2], G).volume() > -1e-10) continue;
					vec3 lI(0, 0, 0);
					FOR(lv, 3) lI += bc_tri[lv] * tr_uvw[lv];


					out_facet = f;
					/*DEBUG*/I = vec3(0, 0, 0); FOR(lv, 3) I += bc_tri[lv] * f.vertex(lv).pos();
					if (!f.opposite().active())
						hex_facet_emb[Volume::Cell(hex, cur_hex).facet(lf)] = f;
					break;
				};

				if (out_facet == -1) {
					Trace::alert("I was not able to cross a tet... maybe due to numerical precision ?");
					break;
				}
				if (!out_facet.opposite().active()) {
					break;
				}
				cur_tet = out_facet.opposite().cell();
				if (tf.ap.mid != 0 && GPTransitionFunction(out_facet, uvw).ap.mid != 0) plop("Double transition");
				tf = tf.apply(GPTransitionFunction(out_facet, uvw));
			}

		}


		if (objective_hex != -1) {// link both cur_hex with objective_hex
			//cur_hex
			FOR(di_cur, 2)FOR(dj_cur, 2)FOR(dk_cur, 2) {
				int lv_cur = di_cur + 2 * dj_cur + 4 * dk_cur;
				vec3 pos_v_cur = -vec3(1, 1, 1) + 2 * vec3(di_cur, dj_cur, dk_cur);
				FOR(di_obj, 2)FOR(dj_obj, 2)FOR(dk_obj, 2) {
					int lv_obj = di_obj + 2 * dj_obj + 4 * dk_obj;
					vec3 pos_v_obj = -vec3(1, 1, 1) + 2 * vec3(di_obj, dj_obj, dk_obj);
					pos_v_obj = tf.ap.get_mat() * pos_v_obj;
					pos_v_obj += 2 * dir;
					if ((pos_v_obj - pos_v_cur).norm2() == 0)
						ds.merge(hex.vert(cur_hex, lv_cur), hex.vert(objective_hex, lv_obj));
				}
			}

			//plrot[ToolBox(pl).add_segment(pos[cur_hex], pos[objective_hex])] = tf.ap.mid;

		}
		else {// it left the volume on I
			//plrot[ToolBox(pl).add_segment(I, pos[cur_hex])] = tf.ap.mid;
		}

	}
	//if (verbose > 1) Drop(pl, plrot).apply_wireframe("links");
	ToolBox(hex).merge_vertices(ds);
	hex.connect();

}



void compute_embedding(Hexahedra& hex, CellFacetAttribute<int>& hex_facet_emb, PointSetEmbedding& emb, TetBoundary& bound) {
	double ave_tri_size = ToolBox(bound.tri).ave_edge_size();
	HexBoundary hexbound(hex);

	FacetAttribute<int> quad_emb(hexbound.quad);

	// transfert quad embedding from CellFacetAttribute(hex) to FacetAttribute(quad)
	// and fill empty quads
	// NOTE: this must be done here to keep HexBoundary local to this function
	{
		for (auto f_quad : hexbound.quad.iter_facets()) {
			auto f_hex = hexbound.hex_facet(f_quad);
			quad_emb[f_quad] = hex_facet_emb[f_hex];
		}
		Drop(hexbound.quad, quad_emb).apply("quademb");


		bool modified;
		do {
			modified = false;
			for (auto h_quad : hexbound.quad.iter_halfedges()) {
				if (quad_emb[h_quad.facet()] == -1 && quad_emb[h_quad.opposite().facet()] != -1) {
					quad_emb[h_quad.facet()] = quad_emb[h_quad.opposite().facet()];
					modified = true;
				}
			}
		} while (modified);
		Drop(hexbound.quad, quad_emb).apply("quadembfilled");
	}



	for (auto v_quad : hexbound.quad.iter_vertices()) {
		if (!v_quad.halfedge().active()) continue;

		//search for feature edge
		std::vector<int> intersected_features;

		for (auto cir : v_quad.iter_halfedges()) {
			Surface::Facet f_tri_src(emb.tri_emb, bound.tri_facet(quad_emb[cir.facet()]));
			Surface::Facet f_tri_dst(emb.tri_emb, bound.tri_facet(quad_emb[cir.opposite().facet()]));

			// search a feature halfedge crossed by a path from f_tri_src to f_tri_dest 
			{// find a path from f_tri_src to f_tri_dst
				auto local_shortest_path_on_facet = [](Triangles& tri, int src, int dest, CornerAttribute<bool>& feature, std::vector<int>& intersected_features,double ave_tri_size) {
					/*TODO: define what is your output*/CornerAttribute<bool> on_path(tri, false);
					std::map<int, double> dist;
					std::map<int, int> prev;
					auto cmp = [&](int i, int j) {// -1 needs special treatment for getting the first eltement
						if (i == -1) return true;
						if (j == -1) return false;
						return dist[i] < dist[j];
						};
					std::set<int, decltype(cmp)> heap(cmp);
					dist[src] = 0;
					prev[src] = -1;
					heap.insert(src);
					while (true) {
						if (heap.empty()) return ;
						auto cur = heap.upper_bound(-1);
						Surface::Facet f(tri, *cur);
						heap.erase(cur);

						if (f == dest) {
							while (f != src) {
								if (feature[prev[f]]) intersected_features.push_back(prev[f]);
								f = Surface::Halfedge(tri, prev[f]).facet();
							}
							return ;
						}

						vec3 G = Triangle3(f).bary_verts();
						for (auto cir : f.iter_halfedges()) {
							auto opp = cir.opposite();
							if (!opp.active()) continue;
							auto opp_f = opp.facet();
							double new_dist = dist[f] + (G + Triangle3(opp_f).bary_verts()).norm();
							if (feature[cir]) new_dist += 50.*ave_tri_size;
							if (dist.find(opp_f) == dist.end() || dist[opp_f] > new_dist) {
								dist[opp_f] = new_dist;
								prev[opp_f] = cir;
								heap.insert(opp_f);
							}
						}
					}
				};
				local_shortest_path_on_facet(emb.tri_emb, f_tri_src, f_tri_dst, emb.feature, intersected_features, ave_tri_size);
			}

		}

		// embed in surface by default
		emb.set_embedding(v_quad, 2, bound.tri_facet(quad_emb[v_quad.halfedge().facet()]));
		
		// embed on feature edge if available 
		if (!intersected_features.empty())
			emb.set_embedding(v_quad, 1, emb.edge_from_halfedge(intersected_features[0])); 

		// embed on sharp vertex is we can find one
		//if(0)
			if (intersected_features.size() > 2) {
			EdgeAttribute<int> multi_grow(emb.pl_emb, 0);
			std::deque<int> queue;
			FOR(i, intersected_features.size()) {
				auto e = emb.edge_from_halfedge(intersected_features[i]);
				multi_grow[e] += 1 << i;
				multi_grow[e.opposite()] += 1 << i;
				queue.push_back(e.opposite());
				queue.push_back(e);
			}
			
			int best_v = -1;
			FOR(iter, 20) {
				if (queue.empty()) break;
				PolyLine::Edge e(emb.pl_emb, queue.front());
				queue.pop_front();
				int valence = 0;
				bool all_touch = false;
				for (auto cir : e.to().iter_edges()) {
					valence++;
					multi_grow[cir] = multi_grow[cir] | multi_grow[e] ;
					if (multi_grow[cir] + 1 == (1 << intersected_features.size()))
						all_touch = true;
				}
				if (valence > 2 && all_touch) {
					if (best_v ==-1) best_v = e.to();
					if (
						(v_quad.pos() - PolyLine::Vertex(emb.pl_emb, best_v).pos()).norm()
						>
						(v_quad.pos() - e.to().pos()).norm())
						best_v = e.to();
						
				}
			}
			if (best_v!=-1) emb.set_embedding(v_quad, 0, best_v);

		}

		vec3 prev_pos = v_quad.pos();
		emb.project_point(v_quad);
		emb.move_toward_point(v_quad, prev_pos);
	}



}

HexCAD& HexCAD::init_hex_from_uvw(TetBoundary& bound, CornerAttribute<bool>& feature, CellCornerAttribute<vec3>& uvw) {
	Triangles& tri = bound.tri;
	Tetrahedra& tet = bound.tet;
	emb.init_from_triangles(tri, &feature);

	double hex_size = ToolBox(tet).ave_edge_size();
	Drop(tet, uvw).apply("uvw");

	// initialize all cubes + their position (geom+map) + the tet that constains them
	Trace::step("init all cubes");
	CellAttribute<int> cube_tet_id(hex);
	CellAttribute<vec3> cube_pos_uvw(hex);
	CellAttribute<vec3> cube_pos(hex);
	DynamicKNN<3> cube_knn(cube_pos.ptr->data);
	generate_cubes(tet, uvw, hex, cube_tet_id, cube_pos_uvw, cube_pos, cube_knn, hex_size);



	CellFacetAttribute<int> hex_facet_emb(hex, -1);
	//PolyLine pl;
	//EdgeAttribute<int> plrot(pl);
	Trace::step("link cubes");
	link_cubes(tet, uvw, hex, cube_tet_id, cube_pos_uvw, cube_pos, cube_knn, hex_facet_emb, hex_size);

	if (verbose > 1) drop_SJ("Connected");

	compute_embedding(hex, hex_facet_emb, emb, bound);
	drop_SJ("Projected");




	emb.show_emb();
	PointAttribute<double> fit_w(hex, 10);
	FOR(iter,4) {
		plop(iter);
	    LeastSquares ls(3*hex.nverts());

		// anisotrop data fitting
		for (auto v : hex.iter_vertices()) {
			auto axes = emb.constrained_direction(v);
			vec3 prev_pos = v.pos();
			emb.project_point(v);
			emb.move_toward_point(v, prev_pos);
			FOR(a, 3-emb.dim[v]) {

				vec3 n = axes[a];
				LinExpr line = - n* v.pos();
				FOR(d, 3) line += n[d] *X(d + 3 * v);
				ls.add_to_energy(fit_w[v] *line);
			}
			v.pos() = prev_pos;
		}


		for (auto h : hex.iter_halfedges()) FOR(d, 3) {
			// preserve epaisseur dans chaque couche 
			ls.add_to_energy(
				X(d + 3 * h.prev().from()) - X(d + 3 * h.prev().to())
				+
				X(d + 3 * h.next().from()) - X(d + 3 * h.next().to())
			);
			// preserve epaisseur des couches successives
			if (!h.facet().opposite().active()) continue;
			ls.add_to_energy(1*(
				X(d + 3 * h.prev().opposite_f().prev().from())
				+ X(d + 3 * h.opposite_c().next().opposite_f().next().to())
				- 2 * X(d + 3 * h.from())
			)
				);

		}
	    ls.solve();
	    for (auto v : hex.iter_vertices()) FOR(d, 3) v.pos()[d] = ls.value(d + 3 * v);
		CellAttribute<double> SJ(hex, 1);
		ToolBox(hex).eval_quality(SJ);
		for (auto v : hex.iter_vertices()) fit_w[v] = 1000;
		for (auto c : hex.iter_cells()) FOR(lv, 8) 
			fit_w[c.vertex(lv)] = 10. * std::max(0.0, std::min(fit_w[c.vertex(lv)],SJ[c]));
		if (verbose > 1) drop_SJ("Optim light");
		Drop(hex,fit_w).apply("fitw");
	}
	


	return*this;
}





