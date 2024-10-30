#include <framework/benjamin_API.h>
#include <fullhex/gp_basic.h>
#include <fullhex/hextract.h>
#include <volume/hex_edit.h>

using namespace UM::Linear;
namespace BenjaminAPI {

	void smooth_vector_on_tet(TetBoundary& bound, FacetAttribute<int>& flag, CellAttribute<vec3>& dir, bool with_tangent_constraint = false) {
		double w_constraint = 100.;
		Triangles& tri = bound.tri;
		Tetrahedra& tet = bound.tet;

		LeastSquares ls(3 * tet.ncells());
		// boundary constraint 
		for (auto f_tri : tri.iter_facets()) {
			vec3 n = Triangle3(f_tri).normal();
			auto c = bound.tet_facet(f_tri).cell();
			if (flag[f_tri] != 0)
				FOR(d, 3) ls.add_to_energy(w_constraint * (X(d + 3 * c) - flag[f_tri] * n[d]));
			else if (with_tangent_constraint) {
				LinExpr line;
				FOR(d, 3) line += n[d] * X(d + 3 * c);
				ls.add_to_energy(w_constraint * line);
			}
		}

		// smoothness
		for (auto f : tet.iter_facets()) {
			if (!f.opposite().active()) continue;
			FOR(d, 3) ls.add_to_energy(X(d + 3 * f.cell()) - X(d + 3 * f.opposite().cell()));
		}

		ls.solve();
		for (auto c : tet.iter_cells()) FOR(d, 3) dir[c][d] = ls.value(d + 3 * c);
	}

	void integrate_vector_on_tet(TetBoundary& bound, FacetAttribute<int>& flag, PointAttribute<double>& scalar, CellAttribute<vec3>& dir, double iso_spacing_coeff = 1) {
		double w_constraint = 100.;
		Triangles& tri = bound.tri;
		Tetrahedra& tet = bound.tet;
		Trace::step("Integrate the objective gradient");
		{
			LeastSquares ls(tet.nverts());
			ls.fix(0, 0);
			// boundary constraint
			for (auto f_tri : tri.iter_facets()) {
				if (flag[f_tri] == 0) continue;
				for (auto h : f_tri.iter_halfedges())
					ls.add_to_energy(w_constraint * (X(h.from()) - X(h.to())));
			}
			for (auto c : tet.iter_cells()) {
				mat<3, 4> grd = Tetrahedron(c).grad_operator();
				LinExpr grad[3];
				FOR(grad_coord, 3) FOR(lv, 4)
					grad[grad_coord] += grd[grad_coord][lv] * X(c.vertex(lv));
				mat3x3 B{ vec3(0, 0, 0), vec3(1, 0, 0),dir[c].normalized() };

				if (std::abs(B[2] * B[1]) > .9) B[1] = vec3(0, 1, 0);
				B[0] = cross(B[2], B[1]).normalized();
				B[1] = cross(B[2], B[0]).normalized();
				FOR(axe, 3) {
					LinExpr line = B[axe][0] * grad[0] + B[axe][1] * grad[1] + B[axe][2] * grad[2];
					if (axe == 2)  line = iso_spacing_coeff * (line - dir[c] * B[axe]);
					ls.add_to_energy(line);
				}
			}
			ls.solve();
			for (auto v : tet.iter_vertices())   scalar[v] = ls.X[v];
		}
	}



	void hextract(Tetrahedra& tet, PointAttribute<vec3>& U, Hexahedra& hex, int nhex_wanted = 10000) {
		DualContour dc(hex);
		Drop(tet, U).apply("U");
		{// rescale
			double volume = 0;
			for (auto c : tet.iter_cells())
				volume += Tetrahedron(U[c.vertex(0)], U[c.vertex(1)], U[c.vertex(2)], U[c.vertex(3)]).volume();
			if (volume < 0) Trace::abort("The volume of the map to extract is negative");
			double scale = std::pow(double(nhex_wanted) / volume, 1. / 3.);
			for (auto v : tet.iter_vertices()) U[v] *= scale;
		}

		// hextract
		CellCornerAttribute<vec3> U_corner(tet);
		for (auto c : tet.iter_corners()) U_corner[c] = U[c.vertex()];
		//Drop(tet, U).apply("U");
		dc.verbose_level(1)
			.init_hex_from_uvw(tet, U_corner)
			//.untangle()
			.drop_SJ("hextracted mesh");

	}







//    _                _                 _                  _____ _____ 
//   | |              (_)               (_)           /\   |  __ \_   _|
//   | |__   ___ _ __  _  __ _ _ __ ___  _ _ __      /  \  | |__) || |  
//   | '_ \ / _ \ '_ \| |/ _` | '_ ` _ \| | '_ \    / /\ \ |  ___/ | |  
//   | |_) |  __/ | | | | (_| | | | | | | | | | |  / ____ \| |    _| |_ 
//   |_.__/ \___|_| |_| |\__,_|_| |_| |_|_|_| |_| /_/    \_\_|   |_____|
//                   _/ |                                               
//                  |__/     



	void polycubify(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, Hexahedra& hex, int nhex_wanted) {
		Trace::SwitchDropInScope drop_switch(false);
		TetBoundary bound(tet);
		PointAttribute<vec3> U(tet);
		FOR(d, 3) {
			FacetAttribute<int> flag(bound.tri);
			for (auto f : tet.iter_facets()) {
				if (tet_flag[f] == -1) continue;
				if (tet_flag[f] % 3 == d)
					flag[bound.tri_facet(f)] = -1 + 2 * (tet_flag[f] / 3);
				else
					flag[bound.tri_facet(f)] = 0;
			}
			Drop(bound.tri, flag)._skip_value(-2).apply("flag");
			CellAttribute<vec3> dir(tet);
			smooth_vector_on_tet(bound, flag, dir);
			PointAttribute<double> scalar(tet);
			integrate_vector_on_tet(bound, flag, scalar, dir);
			Drop(bound.tri, scalar)._skip_value(-2).apply("scalar");
			for (auto v : tet.iter_vertices()) U[v][d] = scalar[v];
		}
		Drop(tet, tet_flag).apply("input");
		hextract(tet, U, hex, nhex_wanted);
	}

	void integrate(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, PointAttribute<vec<3>>& U) {
		Trace::SwitchDropInScope drop_switch(false);
		TetBoundary bound(tet);
		
		FOR(d, 3) {
			FacetAttribute<int> flag(bound.tri);
			for (auto f : tet.iter_facets()) {
				if (tet_flag[f] == -1) continue;
				if (tet_flag[f] % 3 == d)
					flag[bound.tri_facet(f)] = -1 + 2 * (tet_flag[f] / 3);
				else
					flag[bound.tri_facet(f)] = 0;
			}
			Drop(bound.tri, flag)._skip_value(-2).apply("flag");
			CellAttribute<vec3> dir(tet);
			smooth_vector_on_tet(bound, flag, dir);
			PointAttribute<double> scalar(tet);
			integrate_vector_on_tet(bound, flag, scalar, dir);
			Drop(bound.tri, scalar)._skip_value(-2).apply("scalar");
			for (auto v : tet.iter_vertices()) U[v][d] = scalar[v];
		}

	}

	bool pad(Hexahedra& hex, CellFacetAttribute<bool>& pad_face) {
		HexPad padder(hex);
		padder.apply(pad_face); 
		return true;
	};
	
};