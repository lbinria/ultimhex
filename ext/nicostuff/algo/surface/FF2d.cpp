
#include <surface/FF2d.h>




	void FF2D::show_grad(bool with_inverse ) {
		FOR(dir, (with_inverse ? 2 : 1)) {
			FOR(dim, 2) {
				FacetAttribute<vec3> v(m);
				for (auto f : m.iter_facets())
					v[f] = 2. * (double(dir) - .5) * Triangle3(f).grad(vec3(U[f.halfedge(0)][dim], U[f.halfedge(1)][dim], U[f.halfedge(2)][dim]));
				Drop(m, v)._force_length(.005).apply_arrow("ff grad" + std::to_string(dim));
			}
		}
	}

	void FF2D::show_iso(bool with_inverse ) {
		FOR(dir, (with_inverse ? 2 : 1)) {
			FOR(dim, 2) {
				FacetAttribute<vec3> v(m);
				for (auto f : m.iter_facets()) {
					Triangle3 tri  = Triangle3(f);
					vec3 grd = tri.grad(vec3(U[f.halfedge(0)][dim], U[f.halfedge(1)][dim], U[f.halfedge(2)][dim]));

					v[f] = 3. * (double(dir) - .5) * cross(tri.normal(), grd);
				}

				Drop(m, v)._force_radius(.02*ToolBox(m).ave_edge_size()).apply_arrow("ff iso" + std::to_string(dim));
			}
		}
	}

	void FF2D::show_vec2(CornerAttribute<vec2>& uv, FacetAttribute<vec2>& in_vec2, double scale ) {
		FacetAttribute<vec3> v(m);
		scale *= ToolBox(m).ave_edge_size();
		for (auto f : m.iter_facets())
			FOR(dim, 2) v[f] += scale * in_vec2[f][dim] * Triangle3(f).grad(vec3(uv[f.halfedge(0)][dim], uv[f.halfedge(1)][dim], uv[f.halfedge(2)][dim]));
		Drop(m, v).apply_arrow("vec2");
	}

	void FF2D::show_features() { 
		Drop(m, feature)._skip_value(false)._decal(0).apply_half_edge("FF features"); }

	void FF2D::show_singularities() {
		PointAttribute<int> ind(m.points);
		for (auto v : m.iter_vertices())ind[v] = index(v);
		Drop(m.points, ind)._select([&](int v) {return ind[v] != 0; }).apply("FF singu");
	}

	void FF2D::show() {
		show_features();
		show_iso();
		DropSurface(m)._show_vertices(false).apply();
		show_singularities();
	}

	void FF2D::detect_feature_by_dihedral_angle_threshold(double threshold ) {
		for (auto h : m.iter_halfedges()) feature[h] = false;
		for (auto h : m.iter_halfedges()) {
			auto opp = h.opposite();
			if (!opp.active()) feature[h] = true;
			else if (acos(Triangle3(h.facet()).normal() * Triangle3(opp.facet()).normal()) > threshold) {
				feature[h] = true;
				feature[opp] = true;
			}
		}
		Drop(m, feature)._skip_value(false).apply_half_edge();
	}

	void FF2D::estimate_curvature(int nb_smoothing_iters ) {
		std::vector<mat3x3> t(m.nfacets());

		// init tensor per facet
		for (auto h : m.iter_halfedges()) {
			auto opp = h.opposite();
			if (!opp.active()) continue;
			if (feature[h]) continue;
			vec3 n[2] = { Triangle3(h.facet()).normal(),Triangle3(opp.facet()).normal() };
			double angle = acos(n[0] * n[1]);
			vec3 e = h.to().pos() - h.from().pos();
			double s = e.norm() * angle;
			e.normalize();
			mat3x3 M;
			M[0][0] = s * e.x * e.x;
			M[1][1] = s * e.y * e.y;
			M[2][2] = s * e.z * e.z;
			M[0][1] = s * e.x * e.y;
			M[0][2] = s * e.x * e.z;
			M[1][2] = s * e.y * e.z;

			M[1][0] = M[0][1];
			M[2][1] = M[1][2];
			M[2][0] = M[0][2];
			t[h.facet()] += M;
		}

		std::vector<mat3x3> t_accum(m.nfacets());
		for (auto f : m.iter_facets()) {
			std::vector<int> facets;
			std::map<int, int> f2dist;
			f2dist[f] = 0;
			facets.push_back(f);
			FOR(i, facets.size()) {
				t_accum[f] += t[facets[i]];
				Surface::Facet curf(m, facets[i]);
				if (f2dist[curf] > nb_smoothing_iters) continue;
				FOR(lv, 3) {
					auto h = curf.halfedge(lv);
					auto opp = h.opposite();
					if (feature[h] && opp.active()) continue;
					if (f2dist.find(opp.facet()) == f2dist.end()) {
						f2dist[opp.facet()] = f2dist[curf] + 1;
						facets.push_back(opp.facet());
					}
				}
			}
		}

		// get eigen vector & value
		for (auto f : m.iter_facets()) {
			auto [eval, evec] = UM::eigendecompose_symmetric(t_accum[f]);
			FOR(coord, 2) FOR(lv, 3) U[f.halfedge(lv)][coord] = evec.col(coord) * f.vertex(lv).pos();

		}
	}

	void FF2D::init_iso_map(CornerAttribute<vec2>& uv) {
		for (auto f : m.iter_facets()) {
			Triangle3 tri = Triangle3(f);
			auto basis = tri.tangent_basis();
			FOR(lc, 3) uv[f.halfedge(lc)] = (basis * f.halfedge(lc).from().pos()).xy();
		}
		Drop(m, uv).apply_texture("random init FF");
	}

	void FF2D::compute_LS(CornerAttribute<vec2>& metric) {
		CornerAttribute<double> alpha(m);
		for (auto f : m.iter_facets()) {
			Triangle2 tri(metric[f.halfedge(0)], metric[f.halfedge(1)], metric[f.halfedge(2)]);
			vec2 dir = tri.grad(vec3(U[f.halfedge(0)][0], U[f.halfedge(1)][0], U[f.halfedge(2)][0]));
			alpha[f] = atan2(dir[1], dir[0]);
		}

		for (auto h : m.iter_halfedges()) {
			auto opp = h.opposite();
			if (!opp.active()) continue;
			if (opp > h) continue;
			vec2 vh = metric[h.next()] - metric[h];
			vec2 vopp = -(metric[opp.next()] - metric[opp]);
			double delta = std::atan2(mat2x2({ vh,vopp }).det(), vh * vopp);
			LS[h] = std::round((alpha[opp.facet()] - alpha[h.facet()] - delta) / (M_PI / 2.));
			LS[opp] = -LS[h];
		}
		for (auto h : m.iter_halfedges()) LS[h] = (LS[h] + 400000) % 4;
	}

	int FF2D::index(Surface::Vertex v) {
		int ret = 0;
		for (auto h : v.iter_halfedges()) {
			if (!h.opposite().active()) return 0;
			ret += LS[h];
		}
		um_assert(ret > -4000); //;) 
		return (ret + 4000) % 4;
	}
	using namespace UM::Linear;
	void FF2D::extrapolate_FF_from_feature(CornerAttribute<vec2>& metric) {
		LeastSquares solver(2 * m.nfacets());
		for (auto h : m.iter_halfedges()) {
			if (feature[h]) {
				vec2 vh = metric[h.next()] - metric[h];
				double angle = 4. * std::atan2(vh[1], vh[0]);
				vec2 rep(std::cos(angle), std::sin(angle));
				FOR(d, 2) solver.fix(2 * h.facet() + d, rep[d]);
			}
		}
		for (auto h : m.iter_halfedges()) {
			auto opp = h.opposite();
			if (!opp.active()) continue;
			vec2 vh = metric[h.next()] - metric[h];
			vec2 vopp = -(metric[opp.next()] - metric[opp]);
			double delta = 4. * std::atan2(mat2x2({ vh,vopp }).det(), vh * vopp);
			double c = std::cos(delta);
			double s = -std::sin(delta);
			solver.add_to_energy(c*X(2 * h.facet()) + s*X(2 * h.facet() + 1) - X(2 * opp.facet()));
			solver.add_to_energy(-s*X(2 * h.facet()) + c*X(2 * h.facet() + 1) - X(2 * opp.facet() + 1));
		}
		solver.solve();
		for (auto h : m.iter_halfedges()) {
			double angle = .25 * std::atan2(solver.X[2 * h.facet() + 1], solver.X[2 * h.facet()]);
			U[h] = mat2x2({ vec2(std::cos(angle), std::sin(angle)), vec2(-std::sin(angle),std::cos(angle)) }) * metric[h];
		}
		Drop(m, U).apply_texture("extrapolate_FF_from_feature");
	}
	void FF2D::smooth_FF(CornerAttribute<vec2>& metric) {
		LeastSquares ls(2 * m.nfacets()+1);
		auto ls_one = X(2 * m.nfacets());
		for (auto h : m.iter_halfedges()) {
			auto opp = h.opposite();
			if (!opp.active()) continue;
			vec2 vh = metric[h.next()] - metric[h];
			vec2 vopp = -(metric[opp.next()] - metric[opp]);
			double delta = 4. * std::atan2(mat2x2({ vh,vopp }).det(), vh * vopp);
			double c = std::cos(delta);
			double s = -std::sin(delta);
			ls.add_to_energy(c*X(2 * h.facet()) + s*X(2 * h.facet() + 1) - X(2 * opp.facet()));
			ls.add_to_energy(-s*X(2 * h.facet()) + c*X(2 * h.facet() + 1) - X(2 * opp.facet() + 1));
		}

		for (auto f : m.iter_facets()) {
			Triangle2 tri(metric[f.halfedge(0)], metric[f.halfedge(1)], metric[f.halfedge(2)]);
			vec2 dir = tri.grad(vec3(U[f.halfedge(0)][0], U[f.halfedge(1)][0], U[f.halfedge(2)][0]));
			if (dir.norm2() < 1e-10) continue;
			double alpha = 4. * atan2(dir[1], dir[0]);

			ls.add_to_energy(0.01 * (X(2 * f + 0) - std::cos(alpha)*ls_one));
			ls.add_to_energy(0.01 * (X(2 * f + 1) - std::sin(alpha)*ls_one));
		}

		ls.solve();
		for (auto h : m.iter_halfedges()) {
			double angle = .25 * std::atan2(ls.X[2 * h.facet() + 1], ls.X[2 * h.facet()]);
			U[h] = mat2x2({ vec2(std::cos(angle), std::sin(angle)), vec2(-std::sin(angle),std::cos(angle)) }) * metric[h];
		}
	}

	void FF2D::init_map_for_sharp_corners(CornerAttribute<vec2>& iso_metric, CornerAttribute<vec2>& corner_metric) {
		m.compact();

		CornerAttribute<double> sector_angle(m, 0);
		CornerAttribute<double> angle_in_sector(m, 0);
		for (auto h : m.iter_halfedges()) {
			if (!feature[h]) continue;
			auto cir = h;
			do {
				angle_in_sector[cir] = sector_angle[h];
				sector_angle[h] += ToolBox(cir).corner_angle();
				cir = cir.prev().opposite();
			} while (cir.active() && !feature[cir]);
			cir = h;
			do {
				sector_angle[cir] = sector_angle[h];
				cir = cir.prev().opposite();
			} while (cir.active() && !feature[cir]);


		}
		//Drop(m, sector_angle).apply_corner();
		//Drop(m, angle_in_sector).apply_corner();


		// compute stretch
		FacetAttribute<vec2> stretch(m, vec2(0, 0));
		for (auto h : m.iter_halfedges()) {
			if (sector_angle[h] > M_PI / 2.) continue;
			if (sector_angle[h] == 0) continue;
			vec2 start_sector = Geom2::rot(-angle_in_sector[h]) * (iso_metric[h.next()] - iso_metric[h]).normalized();
			vec2 x = Geom2::rot(sector_angle[h] / 2.) * start_sector;
			vec2 y = Geom2::rot(M_PI / 2.) * x;
			double scale = (start_sector * x) / (start_sector * (-y)) - 1.;
			stretch[h.facet()] = scale * x;
		}

		show_vec2(iso_metric, stretch);

		{// smooth stretch 
			LeastSquares solver(2 * m.nfacets());
			for (auto f : m.iter_facets()) if (stretch[f].norm2() != 0)
				FOR(d, 2) solver.fix(2 * f + d, stretch[f][d]);

			for (auto f : m.iter_facets())
				FOR(d, 2) solver.add_to_energy(.1*X(2 * f + d));


			for (auto h : m.iter_halfedges()) {
				if (feature[h]) continue;
				auto opp = h.opposite();
				if (!opp.active()) continue;
				vec2 vh = iso_metric[h.next()] - iso_metric[h];
				vec2 vopp = -(iso_metric[opp.next()] - iso_metric[opp]);
				double delta = std::atan2(mat2x2({ vh,vopp }).det(), vh * vopp);
				double c = std::cos(delta);
				double s = -std::sin(delta);
				solver.add_to_energy(c*X(2 * h.facet()) + s*X(2 * h.facet() + 1) - X(2 * opp.facet()));
				solver.add_to_energy(-s*X(2 * h.facet()) + c*X(2 * h.facet() + 1) - X(2 * opp.facet() + 1));
			}
			solver.solve();
			for (auto f : m.iter_facets()) FOR(d, 2) stretch[f][d] = solver.X[2 * f + d];
		}
		show_vec2(iso_metric, stretch);


		// get new metric
		for (auto f : m.iter_facets()) {
			double scale = stretch[f].norm();
			vec2 x = stretch[f].normalized();
			vec2 y = Geom2::rot(M_PI / 2.) * x;
			if (scale == 0) {
				FOR(lh, 3) corner_metric[f.halfedge(lh)] = iso_metric[f.halfedge(lh)];
				continue;
			}
			FOR(lh, 3) {
				corner_metric[f.halfedge(lh)] = (iso_metric[f.halfedge(lh)] * y) * y + (1. / (scale + 1.)) * (iso_metric[f.halfedge(lh)] * x) * x;
			}
		}
		Drop(m, corner_metric).apply_texture("FF corner_metric");
	}

	void FF2D::use_case_CAD() {
		//Trace::SwitchDropInScope nodrop(false);
		detect_feature_by_dihedral_angle_threshold();
		use_case_CAD_feature_is_given();
	}
	void FF2D::use_case_CAD_feature_is_given() {
		//Trace::SwitchDropInScope nodrop(false);
		CornerAttribute<vec2> iso_metric(m);
		init_iso_map(iso_metric);
		CornerAttribute<vec2> corner_metric(m);
		init_map_for_sharp_corners(iso_metric, corner_metric);
		extrapolate_FF_from_feature(corner_metric);
		compute_LS(corner_metric);
	}
	void FF2D::use_case_scanned_mesh() {
		Trace::SwitchDropInScope nodrop(false);
		estimate_curvature();
		CornerAttribute<vec2> iso_metric(m);
		init_iso_map(iso_metric);
		smooth_FF(iso_metric);
		compute_LS(iso_metric);
	}
