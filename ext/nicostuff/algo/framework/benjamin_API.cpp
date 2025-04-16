#include <framework/benjamin_API.h>
#include <fullhex/gp_basic.h>
#include <fullhex/hextract.h>
#include <volume/hex_edit.h>

#include <surface/feature_curve_detector.h>
#include <surface/pointset_in_surface.h>


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



	void hextract(Tetrahedra& tet, PointAttribute<vec3>& U, Hexahedra& hex, int nhex_wanted, CellFacetAttribute<int> &emb_out) {
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
		// CellFacetAttribute<int> emb_out(hex);
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




	void polycubify(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, Hexahedra& hex, int nhex_wanted, CellFacetAttribute<int> &emb_out) {
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
		hextract(tet, U, hex, nhex_wanted, emb_out);


		Trace::step("Generate tet mesh and transfert flag");
		{
			KNN<3> knn(*(bound.tri.points.data));
			Trace::step("transfert emb from tet facet to tri facet");
			for (auto hex_f : hex.iter_facets()) {
				if (!hex_f.on_boundary()) {
					emb_out[hex_f] = -1;
					continue;
				}
				Volume::Facet tet_f(tet, emb_out[hex_f]);
				int tri_verts[3];
				FOR(lv, 3) tri_verts[lv] = knn.query(tet_f.vertex(lv).pos(), 1)[0];
				for (auto h : Surface::Vertex(bound.tri, tri_verts[0]).iter_halfedges()) {
					if (h.to() == tri_verts[1] && h.prev().from() == tri_verts[2])
						emb_out[hex_f] = h.facet();
				}
			}
		}

		{ // fill empty emb by neigborgs
			HexBoundary hexbound(hex);
			bool done = false;

			while (!done) {
				done = true;
				for (auto h : hexbound.quad.iter_halfedges()) {
					auto opp = h.opposite();
					if (!opp.active()) { Trace::alert("non manifold boundary detected"); continue; }
					if (emb_out[hexbound.hex_facet(h.facet())] == -1) {
						emb_out[hexbound.hex_facet(h.facet())] = emb_out[hexbound.hex_facet(opp.facet())];
						done = false;
					}
				}
			}
		}

	}

	void polycubify(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, Hexahedra& hex, int nhex_wanted) {
		CellFacetAttribute<int> emb_out(hex, -1);
		polycubify(tet, tet_flag, hex, nhex_wanted, emb_out);
	}

	bool pad(Hexahedra& hex, CellFacetAttribute<bool>& pad_face, CellFacetAttribute<int> &emb_attr) {
		HexPad padder(hex);
		padder.apply(pad_face); 

		// Embedding propagate
		{
			for (auto f : hex.iter_facets()) {
				if (!pad_face[f])
					continue;
				
				auto opp = f.halfedge(0).opposite_c();
				if (opp.active()) {
					auto opp_f = opp.opposite_f().next().next().opposite_f().facet();
					std::swap(emb_attr[f], emb_attr[opp_f]);
				}

				for (auto h : f.iter_halfedges()) {
					auto x = h.opposite_f().facet();
					auto opp = h.opposite_c();
					if (opp.active()) {
						auto adj_f = opp.opposite_f().facet();
						
						if (adj_f.on_boundary()) {
							// Copy embedding of adjacent facet
							emb_attr[adj_f] = emb_attr[x];
						}
					}
				}
			}
		}

		return true;
	}

	void embeditinit(Triangles &tri, FacetAttribute<int> &tri_chart, Hexahedra &hex, CellFacetAttribute<int> &emb_attr, FacetAttribute<int> &quad_chart, bool gmsh_chart) {
		// Triangles tri;
		// auto tri_attr = read_by_extension(path + "/tri.geogram", tri);
		tri.connect();

		// Hexahedra hex;
		// auto attr = read_by_extension(path + "/hex.geogram", hex);
		// CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
		hex.connect();

		// FacetAttribute<int> tri_chart(tri);
		// if (gmsh_chart) {
		// 	FacetAttribute<int> gmsh_tri_chart("region", tri_attr, tri);
		// 	for (auto f : tri.iter_facets()) tri_chart[f] = gmsh_tri_chart[f];
		// }
		// else
		// {
			CornerAttribute<int> feature(tri);
			// FeatureEdgeDetector(tri).dihedral_angle().threshold().remove_small_features().remove_small_features().remove_small_features().apply(feature, false);
			// FeatureEdgeDetector(tri).dihedral_angle().threshold().remove_small_features().remove_small_features().remove_small_features().apply(feature, false);
			FeatureEdgeDetector(tri).dihedral_angle().threshold().apply(feature, false);
			// FeatureEdgeDetector(tri).dihedral_angle().threshold().remove_small_features().apply(feature, false);
		

			// Drop(tri, feature)._wireframe(true).apply_half_edge("features");

			DisjointSet ds(tri.nfacets());
			for (auto h : tri.iter_halfedges()) {
				auto opp = h.opposite();
				um_assert(opp.active());
				if (feature[h] == -1)
					ds.merge(h.facet(), opp.facet());
			}
			ds.get_sets_id(tri_chart.ptr->data);
		// }
		// Drop(tri, tri_chart).apply("tri_chart");

		HexBoundary bound(hex);
		Quads& quad = bound.quad;
		// FacetAttribute<int> quad_chart(quad);
		for (auto f : quad.iter_facets()) 
			quad_chart[f] = tri_chart[emb_attr[bound.hex_facet(f)]];
		
		// Drop(quad, quad_chart).apply("tri_chart");

		// if (run_from_graphite) {
		// 	DropSurface(tri).add(tri_chart, "chart")._just_save_filename(path + "/trichart.geogram").apply();
		// 	DropSurface(quad).add(quad_chart, "chart")._just_save_filename(path + "/quadchart.geogram").apply();
		// }
	}

	void embeditapply(Hexahedra &hex, CellFacetAttribute<int> &emb_attr, Quads &quad, FacetAttribute<int> &quad_chart, Triangles &tri, FacetAttribute<int> &tri_chart) {

		FacetAttribute<int> quad2hex_face(quad,-1);
		{
			for (auto f : hex.iter_facets()) {
				if (!f.on_boundary()) continue;
				Surface::Vertex v(quad, f.halfedge(0).from());
				for (auto h : v.iter_halfedges()) if (h.to() == f.halfedge(0).to())
					quad2hex_face[h.facet()] = f;
			}
			for (auto f : quad.iter_facets()) um_assert(quad2hex_face[f] != -1);
		}

		FacetAttribute<int> quad2tri_face(quad, -1);
		FacetAttribute<vec3> bary_pos(quad);
		{
			PolyLine pl;
			for (auto f_quad : quad.iter_facets()) {
				int chart = quad_chart[f_quad];
				vec3 G= Quad3(f_quad).bary_verts();
				quad2tri_face[f_quad] = -1;
				double best_dist2 = 1e20;
				for (auto f_tri : tri.iter_facets()) {
					if (tri_chart[f_tri] != chart) continue;
					vec3 bc = Triangle3(f_tri).bary_coords(G);
					FOR(lv, 3) bc[lv] = std::max(.001, std::min(.999, bc[lv]));
					double sum = bc[0] + bc[1] + bc[2];
					FOR(lv, 3) bc[lv] /= sum;
					vec3 proj = bc[0] * f_tri.vertex(0).pos() + bc[1] * f_tri.vertex(1).pos() + bc[2] * f_tri.vertex(2).pos();
					if ((proj - G).norm2() < best_dist2) {
						quad2tri_face[f_quad] = f_tri;
						bary_pos[f_quad] = proj;
						best_dist2 = (proj - G).norm2();
					}
				}
				ToolBox(pl).add_segment(G, bary_pos[f_quad]);
				um_assert(quad2tri_face[f_quad] != -1);
				emb_attr[quad2hex_face[f_quad]] = quad2tri_face[f_quad];
			}
			DropPolyLine(pl).apply("match");
		}

		// B. comment
		// save_hex_if_valid(hex, emb_attr);

		return;

		// B. comment
		// Drop(quad, quad_chart).apply("quadchart");

		enum {BARY,DIFFUSION,LSCM
		} proj_strat = LSCM;
		switch (proj_strat) {
		case BARY:
			for (auto v : quad.iter_vertices()) {
				vec3 sum_P(0, 0, 0);
				double n = 0;
				for (auto h : v.iter_halfedges()) {
					sum_P += bary_pos[h.facet()];
					n += 1;
				}
				v.pos() = sum_P / n;
			}
			break;
		case (DIFFUSION):
			FOR(dim, 3) {
			double t = 10;
			LeastSquares ls(quad.nverts());
				for (auto v : quad.iter_vertices()) {
					LinExpr line;
					if (!v.halfedge().active()) continue;
					for (auto h : v.iter_halfedges())
						line += X(h.to()) - X(v);
					ls.add_to_energy((X(v) - t * line) - v.pos()[dim]);
				}
				ls.solve();
				for (auto v : quad.iter_vertices()) v.pos()[dim] = ls.value(v);
			}
			break;
		case (LSCM):
			{
				LeastSquares ls(3*quad.nverts());
				for (auto f_quad : quad.iter_facets()) {
					Surface::Facet f_tri(tri, quad2tri_face[f_quad]);
					vec3 G = Triangle3(f_tri).bary_verts();
					vec3 n = Triangle3(f_tri).normal();
					Quaternion q;
					q.v = n * std::sin(M_PI / 4.);
					q.w = std::cos(M_PI / 4.);
					mat3x3 R = q.rotation_matrix();
					for (auto h : f_quad.iter_halfedges()) {
						//lscm
						auto other = h.prev().opposite();
						FOR(dim, 3) {
							LinExpr line = X(dim + 3 * h.to()) - X(dim + 3 * h.from());
							FOR(d, 3) line += R[dim][d] * (X(d + 3 * other.to()) - X(d + 3 * other.from()));
							ls.add_to_energy(line);
						}

						{
							LinExpr line  = -n *G;
							FOR(d, 3) line += n[d] * X(d + 3 * h.from());
							ls.add_to_energy(10.*line);
						}
					}
					//FOR(lv,4) FOR(d,3) ls.add_to_energy(X(d+3*f_quad.vertex(lv))-G[d]);
				}
				ls.solve();
				for (auto v : quad.iter_vertices()) FOR(d,3)v.pos()[d] = ls.value(d+3*v);
			}
			break;

		}
		// B. comment
		// Drop(quad, quad_chart).apply("quadchartsmooth");
	}

	void check_hex_validity(Hexahedra& hex, CellFacetAttribute<int>& emb, std::string msg) {
		EdgeGraph eg(hex);
		for (auto e : eg.iter_edges()) if (!e.opposite().active()) {
			Trace::abort(msg + ": hex mesh has a non manifold edge");
		}
		for (auto f : hex.iter_facets()) if (f.on_boundary() && emb[f] < 0) {
			Trace::abort(msg + " : hex mesh a facet without embedding");
		}
	}

	struct EmbeddedHexSmoother {
		EmbeddedHexSmoother(
			Hexahedra& hex,
			CellFacetAttribute<int>& emb_attr,
			Triangles& tri,
			FacetAttribute<int>& tri_chart) : hex(hex), emb_attr(emb_attr), tri(tri), tri_chart(tri_chart),
			bound(hex), quad(bound.quad), emb(bary) {
			bary.create_points(quad.nfacets());
			for (auto f : quad.iter_facets()) bary[f] = Quad3(f).bary_verts();
			CornerAttribute<bool> feature(tri, false);
			for (auto h : tri.iter_halfedges()) feature[h] = (tri_chart[h.facet()] != tri_chart[h.opposite().facet()]);
			emb.init_from_triangles(tri, &feature);
			for (auto f : quad.iter_facets())  emb.set_embedding(f, 2, emb_attr[bound.hex_facet(f)]);

		}

		Hexahedra& hex;
		CellFacetAttribute<int>& emb_attr;

		Triangles& tri;
		FacetAttribute<int>& tri_chart;
		HexBoundary bound;
		Quads &quad;
		PointSet bary;
		PointSetEmbedding  emb;


		void bary2verts(std::string name) {
			for (auto v : quad.iter_vertices()) {
				vec3 sum_P(0, 0, 0);
				double n = 0;
				for (auto h : v.iter_halfedges()) {
					sum_P += bary[h.facet()];
					n += 1;
				}
				if (n > 0) v.pos() = sum_P / n;
			}
			DropSurface(quad).apply(name);
		};

		void show_boundary_constraints(){// show boundary constraints
			PointAttribute<vec3> n(bary);
			for (auto f : quad.iter_facets())n[f] = Triangle3(Surface::Facet(tri, emb.id[f])).normal();;
			Drop(bary, n).apply("n");
			DropPointSet(bary).apply("bary");
			DropPointSet(emb.pts).apply("emb.pts");
		}

		void lscm() {
			LeastSquares ls(3 * quad.nverts());
			for (auto f_quad : quad.iter_facets()) {
				Surface::Facet f_tri(tri, emb.id[f_quad]);
				vec3 G = Triangle3(f_tri).bary_verts();
				vec3 n = Triangle3(f_tri).normal();
				Quaternion q;
				q.v = n * std::sin(M_PI / 4.);
				q.w = std::cos(M_PI / 4.);
				mat3x3 R = q.rotation_matrix();
				for (auto h : f_quad.iter_halfedges()) {
					auto other = h.prev().opposite();
					FOR(dim, 3) {
						LinExpr line = X(dim + 3 * h.to()) - X(dim + 3 * h.from());
						FOR(d, 3) line += R[dim][d] * (X(d + 3 * other.to()) - X(d + 3 * other.from()));
						ls.add_to_energy(line);
					}
					{
						LinExpr line = -n * G;
						FOR(d, 3) line += n[d] * X(d + 3 * h.from());
						ls.add_to_energy(10. * line);
					}
				}
			}
			ls.solve();
			for (auto v : quad.iter_vertices()) FOR(d, 3)v.pos()[d] = ls.value(d + 3 * v);
		}


		void smooth_inside() {
			LeastSquares ls(3 * quad.nverts());
			// fitting
			for (auto f : quad.iter_facets()) for (auto h : f.iter_halfedges()) 
				FOR(d,3) ls.fix(d + 3 * h.from(),h.from().pos()[d]);
			// LAPLACIEN == no long edges
			double eps = 1;
			for (auto h : hex.iter_halfedges()) FOR(d, 3)
				ls.add_to_energy(eps * (X(d + 3 * h.to()) - X(d + 3 * h.from())));
			ls.solve();
			for (auto v : quad.iter_vertices()) FOR(d, 3) v.pos()[d] = ls.value(d + 3 * v);
			DropVolume(hex).apply("smoothhex");
		}

		void apply() {
			
			
			FOR(it, 3) {
				emb.project();
				PointAttribute<vec3> dest(bary);
				for (auto f : quad.iter_facets()) dest[f] = Quad3(f).bary_verts();
				emb.move_toward(dest);
				show_boundary_constraints();
				lscm();
				DropSurface(quad).apply("lscm");
			}

			smooth_inside();
			DropVolume(hex).apply("smooth_inside");

			//PointAttribute<vec3> dest(bary);
			//FOR(v, bary.size()) dest[v] = bary[v] + vec3(0, 3, 0);
			//bary2verts("quad_init");
			//emb.project();
			//bary2verts("quad_proj");
			//emb.move_toward(dest);
			//bary2verts("quad_walk");
			//emb.project();
			//bary2verts("quad_reproj");
			//show_boundary_constraints();
		}

	};

	void smooth(Hexahedra &hex, CellFacetAttribute<int>&emb_attr, Triangles &tri, FacetAttribute<int> &tri_chart) {

		EmbeddedHexSmoother smoother(hex, emb_attr, tri, tri_chart);
		smoother.apply();
		
		check_hex_validity(hex, emb_attr, "hex mesh validity test FAILED ");
	}

	
};