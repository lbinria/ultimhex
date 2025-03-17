#include <ultimaille/all.h>
#include <framework/trace.h>

#include <framework/nico_framework.h>
#include <framework/benjamin_API.h>
#include "dirty/readonly_mesh_extract_3d.h"
#include <algo/surface/pointset_in_surface.h>


using namespace UM::Linear;






void GP_smooth_vector_on_tet(TetBoundary& bound, FacetAttribute<int>& flag, CellAttribute<vec3>& dir, int dim) {
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
}


void quantize_dim(Tetrahedra& tet, CellFacetAttribute<int>& flag,PointAttribute<vec3>& U,int dim) {
	TetBoundary bound(tet);
	Triangles& tri = bound.tri;
	std::vector<int> chart(tri.nfacets()+1);
	int mass = tri.nfacets();
	DisjointSet ds(mass +1);
	for (auto f : tri.iter_facets()) if (flag[bound.tet_facet(f)] / 2 != dim) ds.merge(mass, f);
	
	for (auto h : tri.iter_halfedges()) {
		auto opp = h.opposite();
		um_assert(opp.active());
		if (flag[bound.tet_facet(h.facet())] / 2 == dim
			&& flag[bound.tet_facet(opp.facet())] / 2 == dim)
			ds.merge(h.facet(), opp.facet());
	}
	int ncharts = ds.get_sets_id(chart);
	//Drop(tri, chart).apply("charts");

	std::vector< std::pair<int, double> > chart_iso(ncharts);
	for (auto h : tri.iter_halfedges()) 
		chart_iso[chart[h.facet()]] = { chart[h.facet()],U[h.from()][dim] };
	chart_iso[chart[mass]] = { chart[mass],1e10 };
	//plop(dim);
	//std::cerr << "AVT------------\n"; for (auto it : chart_iso) std::cerr << it.first << " " << it.second << std::endl;
	std::ranges::sort(chart_iso, [](const std::pair<int, double>& A, const std::pair<int, double>& B) {return A.second < B.second; });
	//std::cerr << "AP------------\n"; for (auto it : chart_iso) std::cerr << it.first << " " << it.second << std::endl;

	double last_value = -10000;
	double decal = 0;
	for (auto& it : chart_iso) {
		it.second = std::floor(.5 + it.second+ decal);
		while (it.second <= last_value) { decal += 1;  it.second += 1; }
	    last_value = it.second;
	}
	//std::cerr << "QUANT------------\n"; for (auto it : chart_iso) std::cerr << it.first << " " << it.second << std::endl;

	std::vector< int> old2new_pos(ncharts);
	FOR(i, ncharts) old2new_pos[chart_iso[i].first] = i;

	for (auto h : tri.iter_halfedges())
		if (flag[bound.tet_facet(h.facet())] / 2 == dim)
			U[h.from()][dim] = chart_iso[old2new_pos[chart[h.facet()]]].second+.5;
}

void GP_integrate_vector_on_tet(Tetrahedra& tet, CellFacetAttribute<int>& flag,
	PointAttribute<vec3>& U, CellAttribute<mat3x3>& B, bool snap = false) {


	auto flag_dim = [&](Volume::Facet f) {
		um_assert(flag[f] >= 0);
		return flag[f] / 2;

		};
	auto flag_sign = [&](Volume::Facet f) {
		um_assert(flag[f] >= 0);
		return (flag[f] % 2) ? -1. : 1.;
		};

	Trace::step("Integrate the objective gradient");
	FOR(dim, 3) {
		if (snap)  quantize_dim(tet, flag, U, dim);
		ConstrainedLeastSquares ls(tet.nverts());
		// boundary constraint
		for (auto f : tet.iter_facets()) {
			if (!f.on_boundary()) continue;
			if (flag_dim(f) != dim) continue;
			if (snap) {
				for (auto h : f.iter_halfedges())
					ls.add_to_constraints(X(h.from()) - U[h.from()][dim]);
			} else 
			for (auto h : f.iter_halfedges()) ls.add_to_constraints(X(h.from()) - X(h.to()));
		}
		for (auto c : tet.iter_cells()) {
			mat<3, 4> grd = Tetrahedron(c).grad_operator();
			LinExpr grad[3];
			FOR(grad_coord, 3) FOR(lv, 4)
				grad[grad_coord] += grd[grad_coord][lv] * X(c.vertex(lv));
			FOR(axe, 3) {
				LinExpr line = B[c][axe][0] * grad[0] + B[c][axe][1] * grad[1] + B[c][axe][2] * grad[2];
				if (axe == dim)  line = (line - B[c][dim] * B[c][axe]);
				ls.add_to_energy(line);
			}
		}
		ls.solve();
		for (auto v : tet.iter_vertices())   U[v][dim] = ls.value(v);
	}
	
}


void polycube_GP(Tetrahedra& tet, CellFacetAttribute<int>& flag, PointAttribute<vec3>& U, int nhex_wanted) {
	//Trace::SwitchDropInScope drop_switch(false);

	auto flag_dim = [&](Volume::Facet f) {
		um_assert(flag[f] >= 0);
		return flag[f] / 2;
		};

	auto flag_sign = [&](Volume::Facet f) {
		um_assert(flag[f] >= 0);
		return (flag[f] % 2) ? -1. : 1.;
		};

	CellAttribute<mat3x3> B(tet, mat3x3::identity());
	FOR(it, 1) {
		Trace::step("FF smoothing iter " + std::to_string(it));
		LeastSquares ls(3 * tet.ncells());
		for (auto f : tet.iter_facets()) {
			if (!f.on_boundary())continue;
			if (flag[f] < 0) continue;
			if (Triangle3(f).unsigned_area() < 1e-10) continue;
			vec3 n = Triangle3(f).normal();

			vec3 rot = -cross(n, flag_sign(f) * B[f.cell()][flag_dim(f)]);
			if (rot.norm2() < 1e-10) continue;

			double rot_n = rot.norm();
			rot.normalize();
			LinExpr line = -rot_n;
			FOR(d, 3) line += rot[d] * X(3 * f.cell() + d);

			ls.add_to_energy(10. * line);
		}
		for (auto f : tet.iter_facets()) {
			auto opp = f.opposite();
			if (!opp.active()) continue;
			FOR(d, 3) ls.add_to_energy(X(d + 3 * f.cell()) - X(d + 3 * opp.cell()));
		}
		ls.solve();
		for (auto c : tet.iter_cells()) {
			Quaternion q;
			FOR(d, 3) q.v[d] = ls.value(d + 3 * c);
			double alpha = q.v.norm();
			if (alpha < 1e-10) continue;// keep identity
			q.v = std::sin(alpha * .5) * q.v.normalized();
			q.w = std::cos(alpha * .5);
			B[c] = B[c] * q.rotation_matrix().transpose();
		}

		//CellAttribute<vec3> vd(tet);
		//FOR(dim, 3) {
		//	for (auto c : tet.iter_cells()) vd[c] = B[c][dim];
		//	Drop(tet, vd).apply_arrow("vec");
		//}
	}
	GP_integrate_vector_on_tet(tet, flag, U, B, false);

	{// rescale
		double volume = 0;
		for (auto c : tet.iter_cells())
			volume += Tetrahedron(U[c.vertex(0)], U[c.vertex(1)], U[c.vertex(2)], U[c.vertex(3)]).volume();
		if (volume < 0) Trace::abort("The volume of the map to extract is negative");
		double scale = std::pow(double(nhex_wanted) / volume, 1. / 3.);
		for (auto v : tet.iter_vertices()) U[v] *= scale;
		for (auto c : tet.iter_cells()) B[c] *= scale;
	}
	GP_integrate_vector_on_tet(tet, flag, U, B, true);



	CellCornerAttribute<vec3> U_corner(tet);
	for (auto c : tet.iter_corners()) U_corner[c] = U[c.vertex()];
	Drop(tet, U_corner).apply("U");
}








void framework_parameters(NicoFramework& fw) {
	fw.add("int", "nbhex", "1000").description("nb hex expected");
	//fw.add("string", "projectpath", "C:/NICO/data/polycubeprojects/B1");
	//fw.add("string", "projectpath", "C:/NICO/data/polycubeprojects/epaule");
	fw.add("string", "projectpath", "C:/NICO/data/mambo-master/mesh/M9.step");
}
void framework_main(NicoFramework& fw) {

	std::string path = fw["projectpath"];
	Triangles tri;
	auto attr = read_by_extension(path + "/flag.geogram", tri);

	tri.connect();
	bool has_flag = false;
	for (auto a : attr.facets) has_flag = has_flag || (a.first.compare("flag") == 0);
	if (!has_flag) {
		std::cerr << "Input mesh does not have a flagging." << std::endl;
		return;
	}
	FacetAttribute<int> flag("flag", attr, tri);

	for (auto f : tri.iter_facets()) {
		um_assert(flag[f] != -1);
		flag[f] /= 5;
	}
	Drop(tri, flag).apply("flag");



	Tetrahedra tet;
	CellFacetAttribute<int> tet_flag(tet, -1);
	Trace::step("Generate tet mesh and transfert flag");
	KNN<3> knn(*(tri.points.data));
	{
		ToolBox<Triangles>(tri).tetgen(tet, true);

		for (auto tet_f : tet.iter_facets()) {
			if (!tet_f.on_boundary()) continue;
			int tri_verts[3];
			FOR(lv, 3) tri_verts[lv] = knn.query(tet_f.vertex(lv).pos(), 1)[0];
			for (auto h : Surface::Vertex(tri, tri_verts[0]).iter_halfedges()) {
				if (h.to() == tri_verts[1] && h.prev().from() == tri_verts[2]) {
					int val = flag[h.facet()];
					tet_flag[tet_f] = val;
				}
			}
		}
		DropVolume(tet).apply("tet");
		Drop(tet, tet_flag).apply("tet_facet");
	}


	PointAttribute<vec3> U(tet);
	{
		//Trace::SwitchDropInScope drop_switch(false);
		polycube_GP(tet, tet_flag, U, fw["nbhex"]);
	}


	Hexahedra hex;
	CellFacetAttribute<int> emb(hex, -1);
	{
		Trace::Section sec("Generate hex mesh and embedding attribute");
		Trace::SwitchDropInScope drop_switch(false);

		Trace::step("Generate");

		CellCornerAttribute<vec3> U_corner(tet);
		for (auto c : tet.iter_corners()) U_corner[c] = U[c.vertex()];
		ReadOnlyMeshExtract3d xtract(tet, U_corner);
		xtract.apply(hex, emb);

		Trace::step("transfert emb from tet facet to tri facet");
		for (auto hex_f : hex.iter_facets()) {
			if (!hex_f.on_boundary()) {
				emb[hex_f] = -1;
				continue;
			}
			Volume::Facet tet_f(tet, emb[hex_f]);
			int tri_verts[3];
			FOR(lv, 3) tri_verts[lv] = knn.query(tet_f.vertex(lv).pos(), 1)[0];
			for (auto h : Surface::Vertex(tri, tri_verts[0]).iter_halfedges()) {
				if (h.to() == tri_verts[1] && h.prev().from() == tri_verts[2])
					emb[hex_f] = h.facet();
			}
		}

	}

	{ // fill empty emb by neigborgs
		HexBoundary bound(hex);
		bool done = false;
		while (!done) {
			done = true;
			for (auto h : bound.quad.iter_halfedges()) {
				auto opp = h.opposite();
				if (!opp.active()) { Trace::alert("non manifold boundary detected"); continue; }
				if (emb[bound.hex_facet(h.facet())] == -1) {
					emb[bound.hex_facet(h.facet())] = emb[bound.hex_facet(opp.facet())];
					done = false;
				}
			}
		}
	}

	Drop(hex, emb).apply("emb");
	//DropVolume(hex).apply("hex");
	DropVolume(hex).add(emb, "emb")._just_save_filename(path + "/hex.geogram").apply();
}
