#include <algo/toolbox.h>
#include <drop_glyph.h>
#include <drop_attribute.h>
#include <toolbox_triangles.h>


	ToolBoxTetrahedra::ToolBoxTetrahedra(Tetrahedra& tets) :tets(tets) {  }
	
	void ToolBoxTetrahedra::peel(double diheral_angle_treshold ) {
		um_assert(tets.connected());
		bool mesh_is_fine = false;
		while (!mesh_is_fine) {
			mesh_is_fine = true;
			std::vector<bool> to_kill(tets.ncells(), false);
			for (auto h : tets.iter_halfedges()) {
				if (!h.facet().on_boundary()) continue;
				if (!h.opposite_f().facet().on_boundary()) continue;
				if (std::abs(ToolBox(h).dihedral_angle()) > diheral_angle_treshold) continue;
				to_kill[h.cell()] = true;
				mesh_is_fine = false;
			}
			tets.delete_cells(to_kill);
			tets.delete_isolated_vertices();
			tets.connect();
		}
	}


	double ToolBoxTetrahedra::ave_edge_size() {
		double sum = 0;
		FOR(c, tets.ncells()) FOR(lv, 4)FOR(lv_inf, lv) sum += (tets.points[tets.vert(c, lv)] - tets.points[tets.vert(c, lv_inf)]).norm();
		return sum / double(6 * tets.ncells());
	}

	void ToolBoxTetrahedra::copy_from(Tetrahedra& other, bool connect ) {
		tets.points.create_points(other.nverts());
		tets.create_cells(other.ncells());
		FOR(v, tets.nverts())tets.points[v] = other.points[v];
		FOR(c, tets.ncells()) FOR(lv, 4) tets.vert(c, lv) = other.vert(c, lv);
		if (connect) tets.connect();
	}


	Triangle3 ToolBoxTetrahedra::facet_geom(int f) {
		return Triangle3(
			tets.points[tets.facet_vert(f / 4, f % 4, 0)],
			tets.points[tets.facet_vert(f / 4, f % 4, 1)],
			tets.points[tets.facet_vert(f / 4, f % 4, 2)]
		);
	}

	double ToolBoxTetrahedra::get_volume() {
		double res = 0;
		for (auto c : tets.iter_cells())
			res += Tetrahedron(c).volume();
		return res;
	}
	void ToolBoxTetrahedra::drop_volume(std::string name) {
		CellAttribute<double> vol(tets, 0);
		for (auto c : tets.iter_cells()) vol[c] = Tetrahedron(c).volume();
		Drop(tets, vol).apply(name);

	}
	void ToolBoxTetrahedra::drop_angles_quality(std::string name, double threshold) {
		CellAttribute<double> val(tets, M_PI/2.);

		for (auto h : tets.iter_halfedges()) {
			double alpha = ToolBox(h).dihedral_angle();
			val[h.cell()] = std::min(val[h.cell()], alpha);
			val[h.cell()] = std::min(val[h.cell()], M_PI-alpha);
		}
		for (auto c : tets.iter_cells()) if (val[c] > threshold) val[c] = -1;
		Drop(tets, val)._skip_value(-1).apply(name);

	}

	void ToolBoxTetrahedra::drop_tangled_cells(std::string name ) {
		CellAttribute<bool> ok(tets, 1);
		for (auto c : tets.iter_cells()) if (Tetrahedron(c).volume() < 0) ok[c] = false;
		Drop(tets, ok).apply(name);
	}


	vec3 ToolBoxTetrahedra::facet_normal(int c, int lf) {
		vec3 P[3];
		FOR(lv, 3) P[lv] = tets.points[tets.facet_vert(c, lf, lv)];
		vec3 n = cross(P[1] - P[0], P[2] - P[0]);
		double norm = n.norm();
		if (norm < 1e-15) return vec3(0, 0, 0); // prevent NANs
		return n / norm;
	}

	void ToolBoxTetrahedra::boundary(Triangles& tri, bool shared_points) {
		//EC3d conn(tets);
		tets.connect();
		if (shared_points) tri.points = tets.points;
		else ToolBox<PointSet>(tri.points).copy_from(tets.points);

		FOR(c, tets.ncells()) FOR(lf, 4) {
			if (tets.conn->oppf[4 * c + lf] != -1) continue;
			int f = tri.create_facets(1);
			FOR(lv, 3) tri.vert(f, lv) = tets.facet_vert(c, lf, lv);
		}
		if (!shared_points)
			tri.delete_isolated_vertices();
	}

	void ToolBoxTetrahedra::read_best_efforts(std::string filename, VolumeAttributes& attribs,bool delete_isolated_vertices ) {
		//void read_tets_with_best_efforts(std::string filename, Tetrahedra& tet, VolumeAttributes& attribs) {
		attribs = read_by_extension(filename, tets);
		if (tets.cells.empty()) {
			Triangles m;
			SurfaceAttributes trash;
			trash = read_by_extension(filename, static_cast<Triangles&>(m));
			m.connect();
			ToolBox<Triangles>(m).tetgen(tets, true);
		}
		if (tets.ncells() == 0) Trace::abort("cannot load file");
		std::cerr << "Tets read from " << filename << "  => #cells= " << tets.ncells() << std::endl;

		//read_tets_with_best_efforts(filename, tets, attribs);
		if (delete_isolated_vertices ) tets.delete_isolated_vertices();
		tets.connect();
	}
	void ToolBoxTetrahedra::read_best_efforts(std::string filename, bool delete_isolated_vertices ) {
		VolumeAttributes attribs;
		read_best_efforts(filename, attribs, delete_isolated_vertices);
	}

	void ToolBoxTetrahedra::check_tet_is_manifold() {
		if (tets.ncells() == 0)
			throw(std::runtime_error("Mesh is empty, I assume that it worst than non-manifold"));

		EdgeGraph eg(tets);

		Trace::step("Find non manifold edges");
		for (auto e : eg.iter_edges())
			for (auto other : e.from().iter_edges())
				if (e.to() == other.to() && e != other) {
					DropVolume(tets).apply("non manifold mesh");
					drop_triangle(Triangle3(e.from().pos(), e.from().pos(), e.to().pos()), "non manifold edge");
					throw(std::runtime_error("non manifold edge"));
				}

		Trace::step("Find non manifold vertices");
		{
			DisjointSet ds(tets.ncells());
			for (auto v : eg.iter_vertices()) {
				std::fill(ds.m_size.begin(), ds.m_size.end(), 1);
				std::iota(ds.m_ids.begin(), ds.m_ids.end(), 0);
				for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
					auto opp_f = cir.facet().opposite();
					if (!opp_f.active()) continue;
					ds.merge(opp_f.cell(), cir.cell());
				}
				int root_id = -1;
				for (auto e : v.iter_edges())for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
					int cell_root = ds.root(cir.cell());
					if (root_id == -1) root_id = cell_root;
					else if (root_id != cell_root) {
						DropVolume(tets).apply("non manifold mesh");
						drop_triangle(Triangle3(v.pos(), v.pos(), v.pos()), "non manifold vertex");
						throw(std::runtime_error("non manifold vertex"));
					}
				}
			}
		}
	}

