#include <framework/benjamin_API.h>
#include <fullhex/gp_basic.h>
#include <fullhex/hextract.h>
#include <volume/hex_edit.h>


#include "dirty/readonly_mesh_extract_3d.h"

using namespace UM::Linear;
namespace BenjaminAPI {

	void smooth_vector_on_tet(TetBoundary& bound, FacetAttribute<int>& flag, CellAttribute<vec3>& dir, int dim) {
		Triangles& tri = bound.tri;
		Tetrahedra& tet = bound.tet;

		LeastSquares ls(3 * tet.ncells());
		for (auto c : tet.iter_cells()) FOR(d, 3) ls.X[d + 3 * c] = d == dim ? 1 : 0;
		// boundary constraint 
		for (auto f_tri : tri.iter_facets()) {
			vec3 n = Triangle3(f_tri).normal();
			auto c = bound.tet_facet(f_tri).cell();
			if (flag[f_tri] == 0) continue;
			FOR(d, 3) ls.fix(d + 3 * c, flag[f_tri] * n[d]);
		}

		// smoothness
		for (auto f : tet.iter_facets()) {
			if (!f.opposite().active()) continue;
			FOR(d, 3) ls.add_to_energy(X(d + 3 * f.cell()) - X(d + 3 * f.opposite().cell()));
		}

		ls.solve();
		for (auto c : tet.iter_cells()) FOR(d, 3) dir[c][d] = ls.value(d + 3 * c);

#if 0
		{// render the field
			PointSet pts;
			PointAttribute<vec3> field(pts);
			pts.create_points(tet.ncells());
			for (auto c : tet.iter_cells()){
				pts[c] = Tetrahedron(c).bary_verts();
				field[c] = dir[c];
			}
			Drop(pts, field).apply("dir");
		}
#endif
	}

	void integrate_vector_on_tet(TetBoundary& bound, FacetAttribute<int>& flag, PointAttribute<double>& scalar, CellAttribute<vec3>& dir,bool snap ) {
		Triangles& tri = bound.tri;
		Tetrahedra& tet = bound.tet;
		Trace::step("Integrate the objective gradient");
		{
			ConstrainedLeastSquares ls(tet.nverts());
			// boundary constraint
			for (auto f_tri : tri.iter_facets()) {
				if (flag[f_tri] == 0) continue;
				if (snap) {
					for (auto h : f_tri.iter_halfedges()) 
						//if (h.from().halfedge() == h)
							//ls.add_to_constraints(X(h.from()) - (.5 + std::floor(.5 + scalar[h.from()])));
							ls.add_to_constraints(X(h.from()) - scalar[h.from()]);
				}
				else for (auto h : f_tri.iter_halfedges()) ls.add_to_constraints(X(h.from()) - X(h.to()));
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
					if (axe == 2)  line = (line - dir[c] * B[axe]);
					ls.add_to_energy(line);
				}
			}
			ls.solve();
			for (auto v : tet.iter_vertices())   scalar[v] = ls.value(v);
		}
	}



	void hextract(Tetrahedra& tet, PointAttribute<vec3>& U, Hexahedra& hex, int nhex_wanted ) {
		Drop(tet, U).apply("U");

		// hextract
		CellCornerAttribute<vec3> U_corner(tet);
		for (auto c : tet.iter_corners()) U_corner[c] = U[c.vertex()];
		DropVolume(tet).add(U_corner,"U")._just_save_filename("C:\\NICO\\tmp\\polyGP.geogram").apply();
		//for (auto c : tet.iter_corners()) U_corner[c] += vec3(.5, .5, .5);

		// show boundary pb

		PointAttribute<vec3> U_danger(tet);
		for (auto v : tet.iter_vertices()) FOR(d, 3) U_danger[v][d] += .5 ;
		for (auto v : tet.iter_vertices()) FOR(d, 3) U_danger[v][d] = std::abs(U[v][d]  - std::floor(.5+U[v][d]));
		Drop(tet, U_danger).apply("U_danger");


		ReadOnlyMeshExtract3d xtract(tet, U_corner);
		CellFacetAttribute<int> emb_out(hex);
		xtract.apply(hex, emb_out);

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
		//Trace::SwitchDropInScope drop_switch(false);
		TetBoundary bound(tet);
		PointAttribute<vec3> U(tet);

		
		// deformation pre-scale&quantization
		FOR(d, 3) {
			FacetAttribute<int> flag(bound.tri);
			for (auto f : tet.iter_facets()) {
				if (tet_flag[f] == -1) continue;
				if (tet_flag[f] % 3 == d)
					flag[bound.tri_facet(f)] = -1 + 2 * (tet_flag[f] / 3);
				else
					flag[bound.tri_facet(f)] = 0;
			}
			//Drop(bound.tri, flag)._skip_value(-2).apply("flag");
			CellAttribute<vec3> dir(tet);
			Trace::step("smooth_vector_on_tet");
			smooth_vector_on_tet(bound, flag, dir,d);
			PointAttribute<double> scalar(tet);
			Trace::step("integrate");
			integrate_vector_on_tet(bound, flag, scalar, dir, false);
			Trace::step("show");
			Drop(bound.tri, scalar)._skip_value(-2).apply("scalar");
			for (auto v : tet.iter_vertices()) U[v][d] = scalar[v];
		}
		{// rescale
			double volume = 0;
			for (auto c : tet.iter_cells())
				volume += Tetrahedron(U[c.vertex(0)], U[c.vertex(1)], U[c.vertex(2)], U[c.vertex(3)]).volume();
			if (volume < 0) Trace::abort("The volume of the map to extract is negative");
			double scale = std::pow(double(nhex_wanted) / volume, 1. / 3.);
			for (auto v : tet.iter_vertices()) U[v] *= scale;
		}

		// deformation and quantization
		FOR(d, 3) {
			FacetAttribute<int> flag(bound.tri);
			for (auto f : tet.iter_facets()) {
				if (tet_flag[f] == -1) continue;
				if (tet_flag[f] % 3 == d)
					flag[bound.tri_facet(f)] = -1 + 2 * (tet_flag[f] / 3);
				else
					flag[bound.tri_facet(f)] = 0;
			}
			CellAttribute<vec3> dir(tet);
			for (auto c : tet.iter_cells()) dir[c] = Tetrahedron(c).grad(vec4{ U[c.vertex(0)][d], U[c.vertex(1)][d], U[c.vertex(2)][d], U[c.vertex(3)][d]});
			
			for (auto c : tet.iter_cells()) if (dir[c].norm2() < 1e-10) dir[c] = vec3(1e-5, 0, 0);

			PointAttribute<double> scalar(tet);
			for (auto v : tet.iter_vertices()) scalar[v] = U[v][d];

			for (auto f : bound.tri.iter_facets()) {
				if (flag[f] ==0) continue;
				for (auto h : f.iter_halfedges())
					scalar[h.from()] = std::floor(scalar[h.from()])+.5;
			}




			Drop(bound.tri, scalar)._skip_value(-2).apply("input scalar");
			Trace::step("integrate");
			integrate_vector_on_tet(bound, flag, scalar, dir,true);
			Trace::step("show");
			Drop(bound.tri, scalar)._skip_value(-2).apply("integrated scalar");
			for (auto v : tet.iter_vertices()) U[v][d] = scalar[v];
		}

	
		for (auto v : tet.iter_vertices()) std::swap(U[v], v.pos());
		DropVolume(tet).apply("tapotte");
		for (auto v : tet.iter_vertices()) std::swap(U[v], v.pos());
		Drop(tet, tet_flag).apply("input");
		hextract(tet, U, hex, nhex_wanted);
	}

	bool pad(Hexahedra& hex, CellFacetAttribute<bool>& pad_face) {
		HexPad padder(hex);
		padder.apply(pad_face); 
		return true;
	};
	
};