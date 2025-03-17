#include<algo/toolbox_triangles.h>
#include<algo/toolbox_pointset.h>
#include <third_party/tetgen/tetgen.h>



	

void ToolBoxTriangles::copy_from(Triangles& other,bool duplicate_vertices) {
	if (duplicate_vertices) ToolBoxPointSet(m.points).copy_from(other.points);
			else m.points = other.points;
	um_assert(other.nfacets() != 0);
	um_assert(m.nfacets() == 0);
	m.create_facets(other.nfacets());
			FOR(f, other.nfacets()) FOR(lv, 3) m.vert(f, lv) = other.vert(f, lv);
		
}



	double ToolBoxTriangles::ave_edge_size() {
		double sum = 0;
		FOR(f, m.nfacets()) FOR(lv, 3) sum += (m.points[m.vert(f, lv)] - m.points[m.vert(f, (lv + 1) % 3)]).norm();
		return sum / double(3 * m.nfacets());
	}

	vec3 ToolBoxTriangles::barycenter(int f) {
		vec3 res(0, 0, 0);
		FOR(lv, m.facet_size(f)) res += m.points[m.vert(f, lv)];
		return res / double(m.facet_size(f));
	}


	void ToolBoxTriangles::split_triangle_4() {
		// add points
		CornerAttribute<int> n_v(m, -1);
		for (auto h : m.iter_halfedges()) {
			if (n_v[h] != -1) continue;
			int v = m.points.create_points(1);
			m.points[v] = .5 * h.to().pos() + .5 * h.from().pos();
			n_v[h] = v;
			auto opp = h.opposite();
			if (opp.active())				n_v[opp] = v;

		}
		// split facets
		int nb_init_facets = m.nfacets(); 
		for (auto f : m.iter_facets()) {
			if (f == nb_init_facets) break;
			FOR(lh, 3)m.conn->create_facet({ f.halfedge(lh).from(),n_v[f.halfedge(lh)], n_v[f.halfedge((lh + 2) % 3)] });
			m.conn->create_facet({ n_v[f.halfedge(0)], n_v[f.halfedge(1)], n_v[f.halfedge(2)]
				});
			m.conn->active[f] = false;
		}
		m.compact();
	}

	void ToolBoxTriangles::tetgen(Tetrahedra& tet, bool with_inner_vertices ) {
		um_assert(m.connected());
		tetgenio in_mesh, out;
		tetgenio::facet* f;
		tetgenio::polygon* p;

		in_mesh.firstnumber = 1;
		in_mesh.numberofpoints = m.nverts();
		in_mesh.pointlist = new REAL[in_mesh.numberofpoints * 3];
		FOR(v, m.nverts()) FOR(d, 3) in_mesh.pointlist[3 * v + d] = m.points[v][d];

		in_mesh.numberoffacets = m.nfacets();
		in_mesh.facetlist = new tetgenio::facet[in_mesh.numberoffacets];
		/**/in_mesh.facetmarkerlist = new int[in_mesh.numberoffacets];
		FOR(f_id, m.nfacets()) {
			f = &in_mesh.facetlist[f_id];
			f->numberofpolygons = 1;
			f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
			f->numberofholes = 0;
			f->holelist = NULL;
			p = &f->polygonlist[0];
			p->numberofvertices = 3;
			p->vertexlist = new int[p->numberofvertices];
			FOR(lv, 3) p->vertexlist[lv] = m.vert(f_id, lv) + 1;
			in_mesh.facetmarkerlist[f_id] = -1;
		}
		tetgenbehavior param;
		double a_val = std::pow(ToolBoxTriangles(m).ave_edge_size(), 3) / 6.;
		char tmp[1024];
		if (with_inner_vertices)
			sprintf(tmp, "YMpq1.1va%f", 20.*a_val);
			//sprintf(tmp, "YMpq1.1va%f", a_val);
		else std::strcpy(tmp, "Mpv");
		param.parse_commandline(tmp);
		tetrahedralize(&param, &in_mesh, &out);

		tet.points.create_points(out.numberofpoints);
		FOR(i, out.numberofpoints)
			FOR(d, 3) tet.points[i][d] = out.pointlist[3 * i + d];

		tet.create_cells(out.numberoftetrahedra);
		FOR(i, out.numberoftetrahedra) FOR(lv, 4){
			tet.vert(i, lv) = out.tetrahedronlist[4 * i + lv] - 1;
		}
		tet.connect();
	}




	void ToolBoxTriangles::read_best_efforts(std::string filename, SurfaceAttributes& attribs) {
		//read_triangles_with_best_efforts(filename, static_cast<Triangles&>(m), attribs);
	
		
		
		attribs = read_by_extension(filename, static_cast<Triangles&>(m));
	if (m.nfacets() == 0) {
		Tetrahedra tets;
		VolumeAttributes trash;
		trash = read_by_extension(filename, tets);
		plop(tets.nverts());
		ToolBox<PointSet>(m.points).copy_from(tets.points);
		//EC3d fec(tets);
		tets.connect();
		for (auto f : tets.iter_facets()) if (tets.conn->oppf[f] == -1) {
			int nf = static_cast<Triangles&>(m).create_facets(1);
			FOR(lv, 3) m.vert(nf, lv) = tets.facet_vert(f / 4, f % 4, lv);
		};
	}

	if (m.nfacets() == 0) Trace::abort("cannot load file");
	m.delete_isolated_vertices();
	std::cerr << "Triangles read from " << filename << "  => #facets = " << m.nfacets() << std::endl;
		
		m.connect();
	}
	void ToolBoxTriangles::read_best_efforts(std::string filename) {
		SurfaceAttributes attribs;
		read_best_efforts(filename, attribs);
	}



	bool ToolBoxTriangles::make_orientable_manifold() {
		//int nb_vertices = m.points.size();

		for (int iter = 0;; iter++) {
			DisjointSet ds(m.nfacets());
			{
				m.connect();
				for (auto h : m.iter_halfedges()) {
					auto opp = h.opposite();
					if (opp.active()) ds.merge(h.facet(), opp.facet());
				}
				m.disconnect();
			}
			FacetAttribute<int> f2chart(m);
			int nbcharts = ds.get_sets_id(f2chart.ptr->data);
			if (nbcharts == 1) {
				m.connect(); 
				return true;
			}
			FOR(f, m.nfacets()) if (f2chart[f] == 0) {
				std::swap(m.vert(f, 1), m.vert(f, 2));
			}
			
			if (iter == 50) {
				m.connect();
				return false;
			}
		}
	}

	int ToolBoxTriangles::add_triangle(Triangle3 tr) {
		int offv = m.points.create_points(3);
		int offf = static_cast<Triangles&>(m).create_facets(1);
		FOR(lv, 3) m.points[offv + lv] = tr[lv];
		FOR(lv, 3) m.vert(offf, lv) = offv + lv;
		return offf;
	}
	void ToolBoxTriangles::add_triangles(Triangles& m1) {
		int offf = static_cast<Triangles&>(m).create_facets(m1.nfacets());
		int offv = m.points.create_points(m1.nverts());
		FOR(f, m1.nfacets()) FOR(lv, 3) m.vert(offf + f, lv) = offv + m1.vert(f, lv);
		FOR(v, m1.nverts()) m.points[offv + v] = m1.points[v];
	}


	//Geom3::Triangle triangle_geom(int f) {
	//	return Geom3::Triangle({ m.points[m.vert(f, 0)], m.points[m.vert(f, 1)], m.points[m.vert(f, 2)] });
	//}



	void ToolBoxTriangles::merge_vertices(double epsilon ) {
		std::vector<int> old2new;
		UM::colocate(*(m.points.data), old2new, epsilon);
		FOR(f, m.nfacets()) FOR(lv, 3) m.vert(f, lv) = old2new[m.vert(f, lv)];
		std::vector<bool> tokill(m.nverts(), true);
		FOR(f, m.nfacets()) FOR(lv, 3) tokill[m.vert(f, lv)] = false;
		m.points.delete_points(tokill, old2new);
		FOR(f, m.nfacets()) FOR(lv, 3) m.vert(f, lv) = old2new[m.vert(f, lv)];
	}

	void ToolBoxTriangles::kill_degenerated_facets_topo() {
		std::vector<bool> to_kill(m.nfacets(), false);
		FOR(f, m.nfacets()) FOR(e, 3) if (m.vert(f, e) == m.vert(f, (e + 1) % 3)) to_kill[f] = true;
		m.delete_facets(to_kill);
	}

	void ToolBoxTriangles::kill_opposite_facets() {
		std::vector<bool> to_kill(m.nfacets(), false);
		FOR(f, m.nfacets()) FOR(e, 3) if (m.vert(f, e) == m.vert(f, (e + 1) % 3)) to_kill[f] = true;

		std::map<std::array<int, 3>, int > map;
		FOR(f, m.nfacets()) {
			std::array<int, 3> key = { m.vert(f, 0),m.vert(f, 1),m.vert(f, 2) };
			std::sort(key.begin(), key.end());
			if (map.end() == map.find(key)) map[key] = f;
			else {
				to_kill[f] = true;
				to_kill[map[key]] = true;
			}
		}
		m.delete_facets(to_kill);
	}



		ToolBoxTriangles::TriangleConnectivityStat::TriangleConnectivityStat(Surface& m) {
			nb_facets = m.nfacets();
			nb_vertices = m.nverts();
			um_assert(m.connected());
			nb_isolated_vertices = 0;
			for (auto v : m.iter_vertices()) if (!v.halfedge().active()) nb_isolated_vertices++;
			for (auto h : m.iter_halfedges()) {
				int nb_opp = 0;
				int nb_dupl = 0;
				for (auto cir : h.from().iter_halfedges()) {
					if (cir.prev().from() == h.to())		nb_opp++;
					if (cir != h && h.to() == cir.to())	nb_dupl++;
				}

				if (nb_opposites.find(nb_opp) == nb_opposites.end()) nb_opposites[nb_opp] = 1;
				nb_opposites[nb_opp]++;
				if (nb_duplicated_edges.find(nb_dupl) == nb_duplicated_edges.end()) nb_duplicated_edges[nb_dupl] = 1;
				nb_duplicated_edges[nb_dupl]++;
			}

			nb_samosas = 0;
			for (auto h : m.iter_halfedges()) {
				auto opp = h.opposite();
				if (!opp.active()) continue;
				if (opp.next().to() == h.next().to())nb_samosas++;
			}
			nb_null_edges = 0;
			for (auto h : m.iter_halfedges()) if (h.to() == h.from()) nb_null_edges++;


			for (auto h : m.iter_halfedges()) for (auto seed : h.from().iter_halfedges()) {
				bool h_is_in_ombrella = false;
				auto cir = seed;
				do {
					h_is_in_ombrella = h_is_in_ombrella || (cir == h);
					cir = cir.prev().opposite();
					if (!cir.active()) {
						h_is_in_ombrella = h_is_in_ombrella || seed.opposite().active();
						break;
					}
				} while (cir != seed);
				if (!h_is_in_ombrella)
					multiple_ombrella_vertices.push_back(h.from());
			}


		}
		std::string ToolBoxTriangles::TriangleConnectivityStat::get_warning_string() {
			std::string result;
			if (nb_vertices == 0)
				result += "#vertices = " + std::to_string(nb_vertices) + "\n";
			if (nb_facets == 0)
				result += "#tri		 = " + std::to_string(nb_facets) + "\n";
			if (nb_isolated_vertices > 0)
				result += "#nb_isolated_vertices		 = " + std::to_string(nb_isolated_vertices) + "\n";

			result += "------------------------------------------------------\n";
			for (auto i : nb_opposites)
				if (i.first != 1)
					result += "We found " + std::to_string(i.second) + " halfedge(s) with " + std::to_string(i.first) + " opposites" + (i.first == 0 ? " (boundary)" : (i.first == 1 ? " (manifold)" : " (non manifold)")) + "\n";
			result += "------------------------------------------------------\n";
			for (auto i : nb_duplicated_edges)
				result += "We found " + std::to_string(i.second) + " halfedge(s) duplicated " + std::to_string(i.first) + " times" + (i.first == 0 ? " (expected)" : " (non manifold)") + "\n";
			result += "nb_samosas = " + std::to_string(nb_samosas) + (nb_samosas == 0 ? " (expected)" : " (WARNING)") + "\n";
			result += "nb_null_edges = " + std::to_string(nb_null_edges) + (nb_null_edges == 0 ? " (expected)" : " (WARNING)") + "\n";
			for (int vid : multiple_ombrella_vertices)
				result += "vertex id = " + std::to_string(vid) + " have mutiple ombrellas (WARNING) \n";
			return result;
		};
		std::string ToolBoxTriangles::TriangleConnectivityStat::get_string() {
			std::string result;
			result += "#vertices = " + std::to_string(nb_vertices) + "\n";
			result += "#tri		 = " + std::to_string(nb_facets) + "\n";
			result += "#nb_isolated_vertices		 = " + std::to_string(nb_isolated_vertices) + "\n";

			result += "------------------------------------------------------\n";
			for (auto i : nb_opposites)
				result += "We found " + std::to_string(i.second) + " halfedge(s) with " + std::to_string(i.first) + " opposites" + (i.first == 0 ? " (boundary)" : (i.first == 1 ? " (manifold)" : " (non manifold)")) + "\n";
			result += "------------------------------------------------------\n";
			for (auto i : nb_duplicated_edges)
				result += "We found " + std::to_string(i.second) + " halfedge(s) duplicated " + std::to_string(i.first) + " times" + (i.first == 0 ? " (expected)" : " (non manifold)") + "\n";
			result += "nb_samosas = " + std::to_string(nb_samosas) + (nb_samosas == 0 ? " (expected)" : " (WARNING)") + "\n";
			result += "nb_null_edges = " + std::to_string(nb_null_edges) + (nb_null_edges == 0 ? " (expected)" : " (WARNING)") + "\n";
			for (int vid : multiple_ombrella_vertices)
				result += "vertex id = " + std::to_string(vid) + " have mutiple ombrellas (WARNING) \n";
			return result;
		};




	void ToolBoxTriangles::laplasmooth(double coeff , bool lock_non_manifold_edges ) {
		um_assert(m.connected());
		std::vector<vec4> tmp(m.nverts()); FOR(v, m.nverts()) FOR(d, 4) tmp[v][d] = 0;
		for (auto h : m.iter_halfedges()) {

			FOR(d, 3) tmp[h.from()][d] += m.points[h.to()][d];
			tmp[h.from()][3]++;
		}
		std::vector<bool> lock(m.nverts(), false);
		if (lock_non_manifold_edges) for (auto h : m.iter_halfedges()) if (!h.opposite().active()) lock[h.from()] = true;

		FOR(v, m.nverts()) if (!lock[v] && tmp[v][3] != 0)
			FOR(d, 3) m.points[v][d] = coeff * m.points[v][d] + (1. - coeff) * tmp[v][d] / tmp[v][3];
	}




