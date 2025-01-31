
#include<algo/toolbox.h>

	ToolBoxHexahedra::ToolBoxHexahedra(Hexahedra& hex) :hex(hex) {  }

	double ToolBoxHexahedra::ave_edge_size() {
		double sum = 0;
		FOR(c, hex.ncells()) FOR(lf, 6) FOR(lv, 4) sum += (hex.points[hex.facet_vert(c, lf, lv)] - hex.points[hex.facet_vert(c, lf, (lv + 1) % 3)]).norm();
		return sum / double(24 * hex.ncells());
	}

	void ToolBoxHexahedra::copy_from(Hexahedra& other, bool duplicate_vertices ) {
		if (duplicate_vertices) ToolBoxPointSet(hex.points).copy_from(other.points);
		else					hex.points = other.points;
		hex.create_cells(other.ncells());
		FOR(c, hex.ncells()) FOR(lv, 8) hex.vert(c, lv) = other.vert(c, lv);
	}


	Quad3 ToolBoxHexahedra::facet_geom(int f) {
		return Quad3(
			hex.points[hex.facet_vert(f / 6, f % 6, 0)],
			hex.points[hex.facet_vert(f / 6, f % 6, 1)],
			hex.points[hex.facet_vert(f / 6, f % 6, 2)],
			hex.points[hex.facet_vert(f / 6, f % 6, 3)]
		);
	}
	Hexahedron ToolBoxHexahedra::cell_geom(int c) {
		return Hexahedron(
			hex.points[hex.vert(c, 0)], hex.points[hex.vert(c, 1)], hex.points[hex.vert(c, 2)], hex.points[hex.vert(c, 3)],
			hex.points[hex.vert(c, 4)], hex.points[hex.vert(c, 5)], hex.points[hex.vert(c, 6)], hex.points[hex.vert(c, 7)]
		);
	}

	void ToolBoxHexahedra::boundary(Quads& quad, bool shared_points) {
		um_assert(hex.connected());
		if (shared_points) quad.points = hex.points;
		else ToolBox<PointSet>(quad.points).copy_from(hex.points);

		FOR(c, hex.ncells()) FOR(lf, 6) {
			if (hex.conn->oppf[6 * c + lf] != -1) continue;
			int f = quad.create_facets(1);
			FOR(lv, 4) quad.vert(f, lv) = hex.facet_vert(c, lf, lv);
		}
		if (!shared_points)
			quad.delete_isolated_vertices();
	}

	void ToolBoxHexahedra::eval_quality(CellAttribute<double>& SJ) {
		hex.connect();
		for (auto c : hex.iter_cells()) {
			SJ[c] = 1;
			for (auto f : c.iter_facets()) {
				for (auto h : f.iter_halfedges()) {
					mat3x3 J = {
						h.to().pos() - h.from().pos(),
						h.opposite_f().next().to().pos() - h.from().pos(),
						h.prev().from().pos() - h.from().pos()
					};
					double div = 1;
					FOR(i, 3) div = div * J.col(i).norm();
					if (div < 1e-15) { SJ[c] = -1.1; continue; }
					SJ[c] = std::min(SJ[c], J.det() / div);
				}
			}
		}
	}

	void ToolBoxHexahedra::drop_scaled_jacobien(std::string name ) {
		CellAttribute<double> SJ(hex);
		eval_quality(SJ);
		Drop(hex, SJ).apply(name);

	}

	void ToolBoxHexahedra::drop_singular_edges(std::string filename, bool show_hex ){
		EdgeGraph eg(hex);
		EdgeAttribute<int> nh(eg, 0);
		EdgeAttribute<bool> boundary(eg, false);
		for (auto h : hex.iter_halfedges()) {
			//if (h % 1000 == 0) std::cerr << int(h) << " / " << hex.nfacets() * 4 << std::endl;
			//auto e = eg.edge_from_halfedge(h);
			PolyLine::Vertex v(eg, h.from());
			PolyLine::Edge e(eg, -1);
			for (auto it : v.iter_edges()) if (it.to() == h.to()) e = it;

			if (!h.facet().opposite().active()) boundary[e] = true;
			nh[e]++;			
		}
		for (auto e : eg.iter_edges()) if (nh[e] == 2 && boundary[e]) nh[e] = 4;
		Drop(eg, nh)._skip_value(4).apply_wireframe(filename);
		if (show_hex) DropVolume(hex).apply(filename + "_hex");

	}

	void ToolBoxHexahedra::drop_boundary(std::string name) {
		CellFacetAttribute<bool> hide(hex);
		for (auto f : hex.iter_facets()) hide[f] = !f.on_boundary();
		Drop(hex, hide).apply("hexboundary");
	}




	void ToolBoxHexahedra::generate_regular_grid(int nb_subdivision) {
		hex.points.create_points(8);
		hex.create_cells(1);
		FOR(lv, 8) hex.points[lv] = vec3(lv % 2, (lv / 2) % 2, lv / 4);
		FOR(lv, 8) hex.vert(0, lv) = lv;
		FOR(i, nb_subdivision)ToolBoxHexahedra(hex).split();
		hex.connect();
	}
	void ToolBoxHexahedra::generate_grid_with_singu_3(int nb_subdivision) {
		hex.points.create_points(14);
		hex.create_cells(3);
		vec3 pts[7] = { vec3(0, 0, 0) ,vec3(.5, 0, 0) ,vec3(1, 0, 0) ,vec3(.2, .5, 0) ,vec3(.5, .4, 0) ,vec3(.8, .5, 0) ,vec3(.5, 1, 0) };
		FOR(i, 7) 	  hex.points[i] = pts[i];
		FOR(lv, 7) hex.points[7 + lv] = hex.points[lv] + vec3(0, 0, .7);

		int topo[3][4] = { {0,1,3,4},{1,2,4,5},{3,4,6,5} };
		FOR(q, 3)FOR(lv, 4) hex.vert(q, lv) = topo[q][lv];
		FOR(c, 3)FOR(lv, 4) hex.vert(c, lv + 4) = hex.vert(c, lv) + 7;
		FOR(i, nb_subdivision)ToolBoxHexahedra(hex).split();
		hex.connect();
	}
	void ToolBoxHexahedra::generate_grid_with_singu_5(int nb_subdivision) {
		hex.points.create_points(22);
		hex.create_cells(5);
		vec3 pts[11] = {
			vec3(.1, .1, 0) ,vec3(.45, .1, 0) ,vec3(.8, 0, 0) ,
			vec3(.1, .45, 0) ,vec3(.5, .5, 0) ,vec3(.9, .3, 0) ,
			vec3(0, .8, 0),vec3(1, .6, 0),
			vec3(.3, .9, 0),vec3(.75, .75, 0),
			vec3(.5, 1, 0)
		};
		FOR(i, 11) 	  hex.points[i] = pts[i];
		FOR(lv, 11) hex.points[11 + lv] = hex.points[lv] + vec3(0, 0, .7);

		int topo[5][4] = { {0,1,3,4},{1,2,4,5},{3,4,6,8},{4,5,9,7},{4,9,8,10} };
		FOR(q, 5)FOR(lv, 4) hex.vert(q, lv) = topo[q][lv];
		FOR(c, 5)FOR(lv, 4) hex.vert(c, lv + 4) = hex.vert(c, lv) + 11;
		FOR(i, nb_subdivision)ToolBoxHexahedra(hex).split();
		hex.connect();
	}




	void ToolBoxHexahedra::split_all(int nb_iter ) {
		if (nb_iter > 1)split_all(nb_iter - 1);
		int off_c = hex.ncells();
		int off_v = hex.points.create_points(64 * hex.ncells());
		hex.create_cells(8 * hex.ncells());
		int c = off_c;
		FOR(old_c, off_c) {
			FOR(off_i, 2)FOR(off_j, 2)FOR(off_k, 2) {
				FOR(di, 2)FOR(dj, 2)FOR(dk, 2) {
					int lv = di + 2 * dj + 4 * dk;
					int v = off_v + 8 * (c - off_c) + lv;
					hex.vert(c, lv) = v;
					double U[3] = { double(off_i + di) / 2.,double(off_j + dj) / 2.,double(off_k + dk) / 2. };
					hex.points[v] = vec3(0, 0, 0);
					FOR(old_di, 2)FOR(old_dj, 2)FOR(old_dk, 2) {
						int old_lv = old_di + 2 * old_dj + 4 * old_dk;
						int old_v = hex.vert(old_c, old_lv);
						double wu = (old_di == 1) ? U[0] : (1. - U[0]);
						double wv = (old_dj == 1) ? U[1] : (1. - U[1]);
						double ww = (old_dk == 1) ? U[2] : (1. - U[2]);
						hex.points[v] += wu * wv * ww * hex.points[old_v];
					}
				}
				c++;
			}
		}
		std::vector<bool> to_kill(off_c, true);
		to_kill.resize(hex.ncells(), false);
		hex.delete_cells(to_kill);

		std::vector<int> old2new;
		ToolBox<PointSet>(hex.points).merge_points(old2new, 1e-10);
		FOR(nc, hex.ncells()) FOR(lv, 8) hex.vert(nc, lv) = old2new[hex.vert(nc, lv)];
		hex.connect();
	}


	/*
	* WARNING: does not support degenerated geometry
	*/
	void ToolBoxHexahedra::split() {
		int off_c = hex.ncells();
		int off_v = hex.points.create_points(64 * hex.ncells());
		hex.create_cells(8 * hex.ncells());
		int c = off_c;
		FOR(old_c, off_c) {
			FOR(off_i, 2)FOR(off_j, 2)FOR(off_k, 2) {
				FOR(di, 2)FOR(dj, 2)FOR(dk, 2) {
					int lv = di + 2 * dj + 4 * dk;
					int v = off_v + 8 * (c - off_c) + lv;
					hex.vert(c, lv) = v;
					double U[3] = { double(off_i + di) / 2.,double(off_j + dj) / 2.,double(off_k + dk) / 2. };
					hex.points[v] = vec3(0, 0, 0);
					FOR(old_di, 2)FOR(old_dj, 2)FOR(old_dk, 2) {
						int old_lv = old_di + 2 * old_dj + 4 * old_dk;
						int old_v = hex.vert(old_c, old_lv);
						double wu = (old_di == 1) ? U[0] : (1. - U[0]);
						double wv = (old_dj == 1) ? U[1] : (1. - U[1]);
						double ww = (old_dk == 1) ? U[2] : (1. - U[2]);
						hex.points[v] += wu * wv * ww * hex.points[old_v];
					}
				}
				c++;
			}
		}
		std::vector<bool> to_kill(off_c, true);
		to_kill.resize(hex.ncells(), false);
		hex.delete_cells(to_kill);

		std::vector<int> old2new;
		ToolBox<PointSet>(hex.points).merge_points(old2new, 1e-10);
		FOR(nvc, hex.ncells()) FOR(lv, 8) hex.vert(nvc, lv) = old2new[hex.vert(nvc, lv)];
		hex.connect();
	}


	void ToolBoxHexahedra::read_best_efforts(std::string filename, VolumeAttributes& attribs) {
		attribs = read_by_extension(filename, hex);
		hex.connect();
	}
	void ToolBoxHexahedra::read_best_efforts(std::string filename) {
		VolumeAttributes attribs;
		read_best_efforts(filename, attribs);
	}

	void ToolBoxHexahedra::merge_vertices(DisjointSet& ds) {
		Trace::step("merge vertices");
		// merge vertices
		PointAttribute<int> set_id(hex.points);
		int new_nverts = ds.get_sets_id(set_id.ptr->data);
		std::vector<vec4> new_pts(new_nverts);
		FOR(v, hex.nverts()) {
			FOR(d, 3) new_pts[set_id[v]][d] += hex.points[v][d];
			new_pts[set_id[v]][3] += 1.;
		}
		FOR(c, hex.ncells()) FOR(lv, 8) hex.vert(c, lv) = set_id[hex.vert(c, lv)];
		FOR(v, new_nverts) hex.points[v] = vec3(new_pts[v][0], new_pts[v][1], new_pts[v][2]) / new_pts[v][3];
		hex.delete_isolated_vertices();
	}



