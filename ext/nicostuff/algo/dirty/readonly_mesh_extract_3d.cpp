#include "dirty/readonly_mesh_extract_3d.h"
#include "toolbox.h"

using namespace Linear;

#define quit(x) {plop(x);Trace::abort();}
inline Tetrahedron uvw_tet(CellCornerAttribute<vec3>& uvw, Volume::Cell c) { return Tetrahedron(uvw[c.corner(0)], uvw[c.corner(1)], uvw[c.corner(2)], uvw[c.corner(3)]); }
inline Triangle3 uvw_tri(CellCornerAttribute<vec3>& uvw, Volume::Facet f) { return Triangle3(uvw[f.corner(0)], uvw[f.corner(1)], uvw[f.corner(2)]); }


vec3 ReadOnlyMeshExtract3d::AxisDirection::dir_vec[6] = { vec3(1,0,0), vec3(0,1,0), vec3(0,0,1),vec3(-1,0,0), vec3(0,-1,0), vec3(0,0,-1) };
int ReadOnlyMeshExtract3d::AxisDirection::dir_inv[6] = { 3,4,5,0,1,2 };


//int transition_rot_times_direction[24][6] = {
//	{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
//	{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
//	{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
//	{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
//	{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
//	{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1}
//};
//FOR(t, 24)FOR(a, 6) FOR(b, 6) {
//	ReadOnlyMeshExtract3d::AxisPermutation dt(t);
//	ReadOnlyMeshExtract3d::AxisDirection da(a);
//	ReadOnlyMeshExtract3d::AxisDirection db(b);
//	if ((dt.matrix() * da.vector() - db.vector()).norm2() == 0)
//		transition_rot_times_direction[t][a] = b;
//}
//FOR(t, 24) {
//	std::cerr << "},\n{";
//	FOR(a, 6)
//		std::cerr << "," << transition_rot_times_direction[t][a];
//}
//return;
int ReadOnlyMeshExtract3d::AxisPermutation::transition_rot_times_direction[24][6] = {
	{ 0, 1, 2, 3, 4, 5 },
	{ 1,3,2,4,0,5 },
{ 3,4,2,0,1,5 },
{ 4,0,2,1,3,5 },
{ 1,2,0,4,5,3 },
{ 3,2,1,0,5,4 },
{ 4,2,3,1,5,0 },
{ 0,2,4,3,5,1 },
{ 2,0,1,5,3,4 },
{ 2,1,3,5,4,0 },
{ 2,3,4,5,0,1 },
{ 2,4,0,5,1,3 },
{ 5,4,3,2,1,0 },
{ 5,0,4,2,3,1 },
{ 5,1,0,2,4,3 },
{ 5,3,1,2,0,4 },
{ 4,3,5,1,0,2 },
{ 0,4,5,3,1,2 },
{ 1,0,5,4,3,2 },
{ 3,1,5,0,4,2 },
{ 3,5,4,0,2,1 },
{ 4,5,0,1,2,3 },
{ 0,5,1,3,2,4 },
{ 1,5,3,4,2,0 } };
//int crossarray[6][6] = { {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1} };
//FOR(a, 6)FOR(b, 6) FOR(c, 6) {
//	ReadOnlyMeshExtract3d::AxisDirection da(a);
//	ReadOnlyMeshExtract3d::AxisDirection db(b);
//	ReadOnlyMeshExtract3d::AxisDirection dc(c);
//	if ((cross(da.vector(), db.vector()) - dc.vector()).norm2() == 0)
//		crossarray[a][b] = c;
//}
//FOR(a, 6) {
//	std::cerr << "},\n{";
//	FOR(b, 6)
//		std::cerr << "," << crossarray[a][b];
//}
//return;
int ReadOnlyMeshExtract3d::AxisDirection::cross[6][6] = {
	{ -1, 2, 4,-1, 5, 1 },
	{  5,-1, 0, 2,-1, 3 },
	{  1, 3,-1, 4, 0,-1 },
	{ -1, 5, 1,-1, 2, 4 },
	{  2,-1, 3, 5,-1, 0 },
	{  4, 0,-1, 1, 3,-1 }
};

int dir2face[6] = { 1,3,5,0,2,4 };
int face2dir[6] = { 3,0,4,1,5,2 };


inline mat<3, 3> mat_from_coeffs(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22) {
	mat<3, 3> res;
	res.rows[0] = vec3(a00, a01, a02);
	res.rows[1] = vec3(a10, a11, a12);
	res.rows[2] = vec3(a20, a21, a22);
	return res;
}

mat3x3 ReadOnlyMeshExtract3d::AxisPermutation::rot_matrix[24] = {
mat_from_coeffs(1,0,0,0,1,0,0,0,1),
	mat_from_coeffs(0,-1,0,1,0,0,0,0,1),mat_from_coeffs(-1,0,0,0,-1,0,0,0,1),mat_from_coeffs(0,1,0,-1,0,0,0,0,1),mat_from_coeffs(0,0,1,1,0,0,0,1,0),
	mat_from_coeffs(-1,0,0,0,0,1,0,1,0),mat_from_coeffs(0,0,-1,-1,0,0,0,1,0),mat_from_coeffs(1,0,0,0,0,-1,0,1,0),mat_from_coeffs(0,1,0,0,0,1,1,0,0),
	mat_from_coeffs(0,0,-1,0,1,0,1,0,0),mat_from_coeffs(0,-1,0,0,0,-1,1,0,0),mat_from_coeffs(0,0,1,0,-1,0,1,0,0),mat_from_coeffs(0,0,-1,0,-1,0,-1,0,0),
	mat_from_coeffs(0,1,0,0,0,-1,-1,0,0),mat_from_coeffs(0,0,1,0,1,0,-1,0,0),mat_from_coeffs(0,-1,0,0,0,1,-1,0,0),mat_from_coeffs(0,-1,0,-1,0,0,0,0,-1),
	mat_from_coeffs(1,0,0,0,-1,0,0,0,-1),mat_from_coeffs(0,1,0,1,0,0,0,0,-1),mat_from_coeffs(-1,0,0,0,1,0,0,0,-1),mat_from_coeffs(-1,0,0,0,0,-1,0,-1,0),
	mat_from_coeffs(0,0,1,-1,0,0,0,-1,0),mat_from_coeffs(1,0,0,0,0,1,0,-1,0),mat_from_coeffs(0,0,-1,1,0,0,0,-1,0)
};
int ReadOnlyMeshExtract3d::AxisPermutation::transition_rot_inverse[24] = {
	0,3,2,1,8,5,15,22,4,14,21,11,12,23,9,6,16,17,18,19,20,10,7,13
};
int ReadOnlyMeshExtract3d::AxisPermutation::transition_rot_times[24][24] = {
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23},
{1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,21,22,23,20},
{2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,18,19,16,17,22,23,20,21},
{3,0,1,2,7,4,5,6,11,8,9,10,15,12,13,14,19,16,17,18,23,20,21,22},
{4,11,21,14,8,3,13,18,0,7,17,22,20,19,5,2,12,23,9,6,16,15,1,10},
{5,8,22,15,9,0,14,19,1,4,18,23,21,16,6,3,13,20,10,7,17,12,2,11},
{6,9,23,12,10,1,15,16,2,5,19,20,22,17,7,0,14,21,11,4,18,13,3,8},
{7,10,20,13,11,2,12,17,3,6,16,21,23,18,4,1,15,22,8,5,19,14,0,9},
{8,22,15,5,0,14,19,9,4,18,23,1,16,6,3,21,20,10,7,13,12,2,11,17},
{9,23,12,6,1,15,16,10,5,19,20,2,17,7,0,22,21,11,4,14,13,3,8,18},
{10,20,13,7,2,12,17,11,6,16,21,3,18,4,1,23,22,8,5,15,14,0,9,19},
{11,21,14,4,3,13,18,8,7,17,22,0,19,5,2,20,23,9,6,12,15,1,10,16},
{12,6,9,23,16,10,1,15,20,2,5,19,0,22,17,7,4,14,21,11,8,18,13,3},
{13,7,10,20,17,11,2,12,21,3,6,16,1,23,18,4,5,15,22,8,9,19,14,0},
{14,4,11,21,18,8,3,13,22,0,7,17,2,20,19,5,6,12,23,9,10,16,15,1},
{15,5,8,22,19,9,0,14,23,1,4,18,3,21,16,6,7,13,20,10,11,17,12,2},
{16,19,18,17,20,23,22,21,12,15,14,13,8,11,10,9,0,3,2,1,4,7,6,5},
{17,16,19,18,21,20,23,22,13,12,15,14,9,8,11,10,1,0,3,2,5,4,7,6},
{18,17,16,19,22,21,20,23,14,13,12,15,10,9,8,11,2,1,0,3,6,5,4,7},
{19,18,17,16,23,22,21,20,15,14,13,12,11,10,9,8,3,2,1,0,7,6,5,4},
{20,13,7,10,12,17,11,2,16,21,3,6,4,1,23,18,8,5,15,22,0,9,19,14},
{21,14,4,11,13,18,8,3,17,22,0,7,5,2,20,19,9,6,12,23,1,10,16,15},
{22,15,5,8,14,19,9,0,18,23,1,4,6,3,21,16,10,7,13,20,2,11,17,12},
{23,12,6,9,15,16,10,1,19,20,2,5,7,0,22,17,11,4,14,21,3,8,18,13}
};
//static int AxisPermutation_stable_direction[24] = {
//2,2,2,2,-1,-1,-1,0,-1,1,-1,-1,-1,-1,1,-1,-1,0,-1,1,-1,-1,0,-1 };


ReadOnlyMeshExtract3d::ReadOnlyMeshExtract3d(Tetrahedra& tet, CellCornerAttribute<vec3>& uvw) :tet(tet), off_hex(tet), uvw(uvw),
link(hex), hex_center_U(hex), tet_cell(hex, -1), emb(hex, { -1,vec3(0,0,0) }), reverted(hex), active(hex, true) {
	double geom_ave_edge_size = ToolBox(tet).ave_edge_size();
	double map_ave_edge_size = 0;
	for (auto h : tet.iter_halfedges()) map_ave_edge_size += (uvw[h.from_corner()] - uvw[h.to_corner()]).norm();
	map_ave_edge_size /= 3 * tet.nfacets();
	ave_edge_size = geom_ave_edge_size / map_ave_edge_size;
}

ReadOnlyMeshExtract3d::AxisDirection cross(ReadOnlyMeshExtract3d::AxisDirection a, ReadOnlyMeshExtract3d::AxisDirection b) {
	return ReadOnlyMeshExtract3d::AxisDirection(ReadOnlyMeshExtract3d::AxisDirection::cross[a][b]);
}
Volume::Facet direction_to_face(Volume::Cell c, ReadOnlyMeshExtract3d::AxisDirection dir) { return c.facet(dir2face[dir]); }
ReadOnlyMeshExtract3d::AxisDirection face_to_direction(Volume::Facet f) { return ReadOnlyMeshExtract3d::AxisDirection(face2dir[f % 6]); }

ReadOnlyMeshExtract3d::AxisDirection halfedge_to_direction(Volume::Halfedge h) { return cross(face_to_direction(h.facet()), face_to_direction(h.opposite_f().facet())); }
Volume::Halfedge direction_to_halfedge(Volume::Facet f, ReadOnlyMeshExtract3d::AxisDirection dir) {
	for (auto h : f.iter_halfedges())
		if (halfedge_to_direction(h).id == dir.id) return h;
	um_assert(!"Should not be there");
}



void ReadOnlyMeshExtract3d::reduce_GP_var(ConstrainedLeastSquares& cls, PointAttribute<bool>& need_jitter) {

	for (auto f : tet.iter_facets()) {
		if (f.on_boundary()) continue;
		if (f.opposite() > f) continue;
		if (!need_jitter[f.vertex(0)] && !need_jitter[f.vertex(1)] && !need_jitter[f.vertex(2)]) continue;
		TransitionFunction tf(f, uvw);
		if (uvw_tri(uvw, f).unsigned_area() < 1e-15) continue;

		tf.force_integer_translation();
		mat3x3 R = tf.m.matrix();
		for (auto h : f.iter_halfedges()) {
			if (!need_jitter[h.from()]) continue;
			Volume::Corner c0 = h.from_corner();
			Volume::Corner c1 = h.opposite_c().to_corner();
			FOR(c0_coord, 3) {
				LinExpr l = (-1) * X(c0_coord + 3 * c0) + tf.t[c0_coord];
				FOR(c1_coord, 3) l = l + R[c0_coord][c1_coord] * X(c1_coord + 3 * c1);	// rotation				
				cls.add_to_constraints(l);
			}
		}
	}
	cls.add_to_energy(LinExpr());
}


void lock_singular_vertices(Tetrahedra& tet, CellCornerAttribute<vec3>& U, PointAttribute<bool>& lock) {
	for (auto v : tet.iter_vertices()) lock[v] = false;
	EdgeGraph eg(tet);
	EdgeAttribute<int> rot(eg, 0);
	EdgeAttribute<bool> on_border(eg, false);
	for (auto e : eg.iter_edges()) {
		auto seed_h = eg.halfedge_from_edge(e);
		ReadOnlyMeshExtract3d::TransitionFunction composed_tf;

		for (auto h : seed_h.iter_CCW_around_edge()) {
			if (!h.facet().opposite().active()) {
				on_border[e] = true;
				continue;
			}
			ReadOnlyMeshExtract3d::TransitionFunction tf(h.facet(), U);
			composed_tf = tf.apply(composed_tf);
		}
		rot[e] = composed_tf.m;
		if (on_border[e]) rot[e] = 0;
		if (rot[e] != 0) {
			lock[e.from()] = true;
			lock[e.to()] = true;
		}
	}
}

void ReadOnlyMeshExtract3d::pre_process() {
	Trace::Section sec("preprocess");

	auto  grid_point_on_triangle = [&](Volume::Facet f) {
		Triangle3 tri(uvw[f.corner(0)], uvw[f.corner(1)], uvw[f.corner(2)]);
		//Triangle3 tri(uvw[f.corner(0)] + vec3(.5, .5, .5), uvw[f.corner(1)] + vec3(.5, .5, .5), uvw[f.corner(2)]+vec3(.5, .5, .5));

		BBox3 box;
		FOR(lv, 3) box.add(tri[lv]);
		box.dilate(1e-2);
		for (int i = std::floor(box.min[0]); i < std::ceil(box.max[0]); i++) {
			for (int j = std::floor(box.min[1]); j < std::ceil(box.max[1]); j++) {
				for (int k = std::floor(box.min[2]); k < std::ceil(box.max[2]); k++) {
					vec3 P = vec3(i, j, k);
					double eps = 1e-5;
					double eps2 = eps * eps;
					// near extremities
					FOR(lv, 3) if ((P - tri[lv]).norm2() < eps2) return true;

					// near edge
					FOR(lh, 3) {
						vec3 A = tri[lh];
						vec3 B = tri[(lh + 1) % 3];
						// not in the segment
						if ((P - A) * (B - A) < 0) continue;
						if ((P - B) * (A - B) < 0) continue;
						if ((A - B).norm2() < eps2) continue;// degenerated edge
						if (cross(B - A, P - A).norm2() / (B - A).norm2() < eps2) return true;
					}
					// near triangle
					vec3 cp = cross(tri[2] - tri[0], tri[1] - tri[0]);
					if (cp.norm2() < eps2) continue;// degenerated triangle
					if (std::abs(cp.normalized() * (P - tri[0])) < eps) {
						vec3 bc = tri.bary_coords(P);
						if (bc[0] < 0) continue;
						if (bc[1] < 0) continue;
						if (bc[2] < 0) continue;
						//if ([&]() {FOR(lv,3) if (bc[lv] < 0) return true; return false; }()) continue;
						return true;
					}
				}
			}
		}
		return false;
		};
	auto  grid_point_on_halfedge = [&](Volume::Halfedge h) {
		Segment3 seg(uvw[h.from_corner()], uvw[h.to_corner()]);
		//Segment3 seg(uvw[h.from_corner()]+vec3(.5,.5,.5), uvw[h.to_corner()] + vec3(.5, .5, .5));
		BBox3 box;
		FOR(lv, 2) box.add(seg[lv]);
		box.dilate(1e-2);
		for (int i = std::floor(box.min[0]); i < std::ceil(box.max[0]); i++) {
			for (int j = std::floor(box.min[1]); j < std::ceil(box.max[1]); j++) {
				for (int k = std::floor(box.min[2]); k < std::ceil(box.max[2]); k++) {
					vec3 P = vec3(i, j, k);
					double eps = 1e-5;
					double eps2 = eps * eps;
					// near extremities
					FOR(lv, 2) if ((P - seg[lv]).norm2() < eps2) return true;

					// near edge
					FOR(lh, 3) {
						vec3 A = seg[0];
						vec3 B = seg[1];
						// not in the segment
						if ((P - A) * (B - A) < 0) continue;
						if ((P - B) * (A - B) < 0) continue;
						if ((A - B).norm2() < eps2) continue;// degenerated edge
						if (cross(B - A, P - A).norm2() / (B - A).norm2() < eps2) return true;
					}
				}
			}
		}
		return false;
		};

	auto show_bad_triangles = [&]() {
		Triangles tri;
		for (auto f : tet.iter_facets()) if (grid_point_on_triangle(f))
			ToolBox(tri).add_triangle({ f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos() });
		DropSurface(tri).apply("grid_point_on_triangle");
		Drop(tet, uvw).apply_neg_det("negdej");
		};







	show_bad_triangles();
	PointAttribute<bool> lock(tet);
	lock_singular_vertices(tet, uvw, lock);
	int last_nfail = std::numeric_limits<int>::max();
	FOR(it, 4) {
		std::cerr << "Jitter iteration " << it << std::endl;
		int nfail = 0;
		PointAttribute<bool> need_jitter(tet, false);
		for (auto f : tet.iter_facets())
			if (grid_point_on_triangle(f))
				for (auto h : f.iter_halfedges()) {
					if (!lock[h.from()]) {
						need_jitter[h.from()] = true;
						nfail++;
					}
				}
		plop(nfail);
		if (nfail == 0 || nfail >= last_nfail) break;
		last_nfail = nfail;
		ConstrainedLeastSquares cls(3 * tet.ncorners());
		reduce_GP_var(cls, need_jitter);
		CRSMatrix& M = cls.M;
		int const_index = M.count_columns() - 1;
		std::vector<double> reduced(M.count_columns(), 1);


		auto change_uvw = [&](int c, int d, double val) {
			auto row = M.row(d + 3 * c);
			if (row.size() != 1) return false;
			if (row[0].index == const_index) return false;
			reduced[row[0].index] = val / row[0].value;
			return true;
			};

		auto propagate_uvw_changes = [&]() {
			for (auto c : tet.iter_corners()) {
				FOR(d, 3) {
					uvw[c][d] = 0;
					for (auto it : M.iter_row(d + 3 * c))
						uvw[c][d] += it.value * reduced[it.index];
				}
			}
			};

		for (auto c : tet.iter_corners()) FOR(d, 3) change_uvw(c, d, uvw[c][d]);
		propagate_uvw_changes();


		for (auto c : tet.iter_corners()) if (need_jitter[c.vertex()])
			FOR(d, 3) change_uvw(c, d, uvw[c][d] + .1 * rand_range(-1, 1));
		//for (auto f : tet.iter_facets()) if (grid_point_on_triangle(f)) for (auto h : f.iter_halfedges()) {
		//	if (lock[h.from()]) continue;
		//	FOR(d, 3) change_uvw(h.from_corner(), d, uvw[h.from_corner()][d] + .1 * rand_range(-1, 1));
		//}
		propagate_uvw_changes();
	}
	show_bad_triangles();


}

Volume::Halfedge ReadOnlyMeshExtract3d::link_opposite(Volume::Halfedge h) {
	Volume::Facet f = h.facet();
	if (!link[f].active()) return  Volume::Halfedge(hex, -1);
	AxisDirection dir_f = face_to_direction(f);
	dir_f = link[f].m * dir_f;

	AxisDirection dir_h = halfedge_to_direction(h);
	dir_h = link[f].m * dir_h;

	if (reverted[f.cell()] == reverted[link[f].q_id]) {
		dir_f = dir_f.inverse();
	}
	else dir_h = dir_h.inverse();

	Volume::Facet fopp = direction_to_face(Volume::Cell(hex, link[f].q_id), dir_f);
	for (auto hopp : fopp.iter_halfedges()) {
		if (halfedge_to_direction(hopp) == dir_h.inverse())
			return hopp;
	}
	um_assert(!"Should be be there");
}

Volume::Facet ReadOnlyMeshExtract3d::link_opposite_facet(Volume::Facet f) {
	if (!link[f].active()) return  Volume::Facet(hex, -1);
	AxisDirection dir_f = face_to_direction(f);
	dir_f = link[f].m * dir_f;

	if (reverted[f.cell()] == reverted[link[f].q_id])  dir_f = dir_f.inverse();

	return direction_to_face(Volume::Cell(hex, link[f].q_id), dir_f);
}

void ReadOnlyMeshExtract3d::smooth_hex_network() {
	CellAttribute<vec3> G(hex);
	for (auto c : hex.iter_cells()) G[c] = Hexahedron(c).bary_verts();
	FOR(d, 3) {
		LeastSquares ls(hex.ncells());
		for (auto f : hex.iter_facets()) if (!link[f].active())
			ls.fix(f.cell(), G[f.cell()][d]);
		for (auto f : hex.iter_facets()) if (link[f].active())
			ls.add_to_energy(X(f.cell()) - X(link[f].q_id));
		ls.solve();
		for (auto c : hex.iter_cells()) FOR(lv, 8)
			c.vertex(lv).pos()[d] += ls.value(c) - G[c][d];
	}
}

void ReadOnlyMeshExtract3d::show_hex_network(bool show_hex, bool show_ribbons) {
	if (show_hex) {
		CellAttribute<int> orient(hex, 0);
		for (auto c : hex.iter_cells()) if (active[c]) orient[c] = reverted[c] ? -1 : 1;  else orient[c] = -2;
		Drop(hex, orient)._skip_value(-2).apply("hex orient");
	}
	//{

	if (show_ribbons) {
		Polygons poly;
		FacetAttribute<int> polyattr(poly, 0);
		double c0 = .5;// shrink coefficient
		double c1 = 1. - c0;
		for (auto h : hex.iter_halfedges()) {
			if (!active[h.cell()]) continue;
			auto opp = link_opposite(h);
			if (!opp.active()) continue;
			vec3 G = Quad3(h.facet()).bary_verts();
			vec3 Gopp = Quad3(opp.facet()).bary_verts();

			if (reverted[h.cell()] == reverted[opp.cell()])
				Glyphs::Ribbon(poly,
					c1 * h.from().pos() + c0 * G, c1 * h.to().pos() + c0 * G,
					G - Hexahedron(h.cell()).bary_verts(),
					c1 * opp.to().pos() + c0 * Gopp, c1 * opp.from().pos() + c0 * Gopp,
					Gopp - Hexahedron(opp.cell()).bary_verts(),
					10).portion_only(0, .5).apply(polyattr, (h.opposite_f().facet()) % 6);
			else
				Glyphs::Ribbon(poly,
					c1 * h.from().pos() + c0 * G, c1 * h.to().pos() + c0 * G,
					G - Hexahedron(h.cell()).bary_verts(),
					c1 * opp.from().pos() + c0 * Gopp, c1 * opp.to().pos() + c0 * Gopp,
					Gopp - Hexahedron(opp.cell()).bary_verts(),
					10).portion_only(0, .5).apply(polyattr, (h.opposite_f().facet()) % 6);
		}
		Drop(poly, polyattr).apply("hex links");
	}

}





using namespace Linear;


Tetrahedron ReadOnlyMeshExtract3d::uvw_tet(Volume::Cell t) { return Tetrahedron(uvw[t.corner(0)], uvw[t.corner(1)], uvw[t.corner(2)], uvw[t.corner(3)]); }

void ReadOnlyMeshExtract3d::create_hex(Volume::Cell t, int i, int j, int k) {
	Tetrahedron tet_U = uvw_tet(t);
	vec4 bc = tet_U.bary_coords(vec3(i, j, k));
	um_assert(bc[0] >= 0 && bc[1] >= 0 && bc[2] >= 0 && bc[3] >= 0);
	int off_f = hex.create_cells(1);
	int off_v = hex.points.create_points(8);
	FOR(lv, 8) hex.vert(off_f, lv) = off_v + lv;
	tet_cell[off_f] = t;
	hex_center_U[off_f] = vec3(i, j, k);
	reverted[off_f] = (uvw_tet(t).volume() < 0);

	// rendering purpose
	{
		Tetrahedron tet_X = Tetrahedron(t);
		vec3 G = (bc[0] * tet_X[0] + bc[1] * tet_X[1] + bc[2] * tet_X[2] + bc[3] * tet_X[3]);
		double s = .5 * ave_edge_size;
		mat3x3 grad = {
			tet_X.grad(vec4{uvw[t.corner(0)][0], uvw[t.corner(1)][0], uvw[t.corner(2)][0], uvw[t.corner(3)][0]}),
			tet_X.grad(vec4{uvw[t.corner(0)][1], uvw[t.corner(1)][1], uvw[t.corner(2)][1], uvw[t.corner(3)][1]}),
			tet_X.grad(vec4{uvw[t.corner(0)][2], uvw[t.corner(1)][2], uvw[t.corner(2)][2], uvw[t.corner(3)][2]})
		};
		grad = grad.invert_transpose();
		FOR(d, 3) grad[d].normalize();
		FOR(di, 2)FOR(dj, 2)FOR(dk, 2)
			hex.points[off_v + dk * 4 + dj * 2 + di] = G + s * (double(di) - .5) * grad[0] + s * (double(dj) - .5) * grad[1] + s * (double(dk) - .5) * grad[2];
	}
}

void ReadOnlyMeshExtract3d::generate_hexes() {
	for (auto t : tet.iter_cells()) {
		if (t % 1000 == 0) std::cerr << "tet " << t << " / " << tet.ncells() << std::endl;
		off_hex[t] = hex.ncells();
		Tetrahedron tet_U = uvw_tet(t);
		if (std::abs(tet_U.volume()) < 1e-10) {
			//plop(tet_U.volume()); 
			continue;
		}
		BBox3 box;
		FOR(lv, 4) box.add(uvw[t.corner(lv)]);
		for (int i = std::floor(box.min[0]); i < std::ceil(box.max[0]); i++) {
			for (int j = std::floor(box.min[1]); j < std::ceil(box.max[1]); j++) {
				for (int k = std::floor(box.min[2]); k < std::ceil(box.max[2]); k++) {
					vec4 bc = tet_U.bary_coords(vec3(i, j, k));
					if (bc[0] < 0 || bc[1] < 0 || bc[2] < 0 || bc[3] < 0) continue;
					//if (std::abs(tet_U.volume()) < 1e-18) continue;//Trace::alert("Hex create on low volume");
					create_hex(t, i, j, k);
				}
			}
		}
	}
	hex.connect();
}


Volume::Cell ReadOnlyMeshExtract3d::search_hex(Volume::Cell t, vec3 ijk) {
	for (int obj_hex = off_hex[t];
		obj_hex < ((t + 1 == tet.ncells()) ? hex.ncells() : off_hex[t + 1]);
		obj_hex++) {
		if ((hex_center_U[obj_hex] - ijk).norm2() < 1e-5)
			return Volume::Cell(hex, obj_hex);
	}
	plop(uvw_tet(t).volume());
	for (int obj_hex = off_hex[t];
		obj_hex < ((t + 1 == tet.ncells()) ? hex.ncells() : off_hex[t + 1]);
		obj_hex++) {
		plop(hex_center_U[obj_hex]);
		plop(ijk);
		plop(uvw_tet(t).bary_coords(ijk));
	}
	um_assert(!"position reached");
}

void ReadOnlyMeshExtract3d::connect_hexes() {
	bool verbose = false;
	for (auto f : hex.iter_facets()) {
		if (verbose) std::cerr << "#######################################################################\n";;
		vec3 out_position(0, 0, 0);
		Volume::Cell seed_tet(tet, tet_cell[f.cell()]);
		vec3 O = hex_center_U[f.cell()];
		Volume::Cell cur_tet = seed_tet;

		TransitionFunction tf;
		Volume::Facet in_facet(tet, -1);
		int niter = 0;
		while (true) {
			bool seed_is_reverted = (uvw_tet(seed_tet).volume() < 0);
			bool cur_is_reverted = (uvw_tet(cur_tet).volume() < 0);
			Tetrahedron tet_U = uvw_tet(cur_tet);
			FOR(lv, 4) tet_U[lv] = tf.apply(tet_U[lv]);

			//if (niter> 99)verbose = true;
			if (niter++ > 100) {
				break;
				plop(Tetrahedron(seed_tet).bary_coords(O));
				plop(cur_tet);
				plop(seed_tet);
				plop(tet_U.volume());
				FOR(lv, 4) plop(tet_U[lv]);
				FOR(lf, 4)drop_triangle(Triangle3(seed_tet.facet(lf)), "failtet");
				drop_point(Hexahedron(f.cell()).bary_verts(), "bary");
				drop_point(out_position, "out_position");

				quit(niter);
				//Trace::abort();
			}


			{// test if  destination is reached

				vec3 dest = O;
				if (seed_is_reverted == cur_is_reverted)
					dest += AxisDirection(face2dir[f % 6]).vector();

				vec4 bc = tet_U.bary_coords(dest);
				if (bc[0] > 0 && bc[1] > 0 && bc[2] > 0 && bc[3] > 0) {
					auto obj_q = search_hex(cur_tet, tf.inverted().apply(dest));
					link[f] = Link(obj_q, tf.m.inverse().id);
					um_assert(link[f].m != -1);
					break;
				}
			}


			vec3 dir = AxisDirection(face2dir[f % 6]).vector();
			if (seed_is_reverted == cur_is_reverted) dir *= -1;

			// cross tet cur_tet
			if (verbose) std::cerr << "iter = " << niter << "\tvol= " << tet_U.volume() << std::endl;

			if (verbose) std::cerr << cur_tet << "/" << seed_tet << std::endl;

			if (verbose)  plop(O);
			if (verbose)  plop(dir);

			bool search_fail = false;
			for (int lf = 0;; lf++) {
				if (lf == 4) {
					break;
					plop(Tetrahedron(seed_tet).bary_coords(O));
					plop(int(f));
					plop(cur_tet);
					plop(seed_tet);
					plop(tet_U.volume());
					FOR(lv, 4) plop(tet_U[lv]);
					FOR(lf, 4)drop_triangle(Triangle3(seed_tet.facet(lf)), "failtet");
					drop_point(Hexahedron(f.cell()).bary_verts(), "bary");
					drop_point(out_position, "out_position");
					//quit(niter);
					search_fail = true;
				};
				//for (auto cir_facet : cur_tet.iter_facets()) {
				Volume::Facet cir_facet = cur_tet.facet(lf);
				Triangle3 tri_U(uvw[cir_facet.corner(0)], uvw[cir_facet.corner(1)], uvw[cir_facet.corner(2)]);
				FOR(lv, 3) tri_U[lv] = tf.apply(tri_U[lv]);

				vec3 w;
				FOR(i, 3) w[i] = Tetrahedron(O, O + dir, tri_U[(i + 1) % 3], tri_U[(i + 2) % 3]).volume();
				if (!cur_is_reverted)FOR(i, 3)w[i] *= -1;

				if (verbose)  plop(w);
				// do not go back
				if (cir_facet == in_facet) continue;

				// path must cross the triangle
				if (w[0] < 0) continue;
				if (w[1] < 0) continue;
				if (w[2] < 0) continue;


				vec3 I = (w[0] * tri_U[0] + w[1] * tri_U[1] + w[2] * tri_U[2]) / (w[0] + w[1] + w[2]);
				if (verbose) plop(w);
				if (verbose) plop(I);

				out_position = (w[0] * cir_facet.vertex(0).pos() + w[1] * cir_facet.vertex(1).pos() + w[2] * cir_facet.vertex(2).pos()) / (w[0] + w[1] + w[2]);
				// goto next facet or reach boundary
				if (cir_facet.opposite().active()) {
					in_facet = cir_facet.opposite();
					TransitionFunction ntf(cir_facet, uvw);
					ntf.force_integer_translation();
					//if (tf.m != 0 && ntf.m != 0) plop(tf.m);
					tf = tf.apply(ntf);
					cur_tet = cir_facet.opposite().cell();
				}
				else {
					emb[f] = { cir_facet,out_position };
				}
				break;
			}
			// boundary is reached
			if (search_fail || emb[f].f_id != -1)break;
		}
	}
	// enforce link symmetry
	for (auto f : hex.iter_facets()) if (link[f].active()) {
		auto opp = link_opposite_facet(f);
		if (link[opp].active() && link_opposite_facet(opp) == f) continue;
		link[f] = Link();
		um_assert(!link[f].active());
	}
	//check that links are symmetric
	for (auto f : hex.iter_facets())
		if (link[f].active()) um_assert(link_opposite_facet(link_opposite_facet(f)) == f);
}



void ReadOnlyMeshExtract3d::propagate_boundary() {
	for (auto f_seed : hex.iter_facets()) {
		if (link[f_seed].active()) continue;
		auto f = f_seed;
		while (f.active()) {

			AxisDirection dir_f = face_to_direction(f);
			f = 2 * (f / 2) + 1 - (f % 2);
			if (!link[f].active()) break;
			dir_f = link[f].m * dir_f;
			if (reverted[f.cell()] != reverted[link[f].q_id])
				dir_f = dir_f.inverse();
			f = direction_to_face(Volume::Cell(hex, link[f].q_id), dir_f);
			um_assert(link[f].q_id != -1);
			um_assert(emb[f].f_id == -1);
			emb[f] = emb[f_seed];
		}
	}
}


bool ReadOnlyMeshExtract3d::can_cancel(Volume::Facet seed_f) {
	auto lnk = link[seed_f];
	if (!lnk.active()) return false;
	Volume::Cell c0 = seed_f.cell();
	Volume::Cell c1(hex, lnk.q_id);
	um_assert(c0 != c1);
	um_assert(active[c0]);
	um_assert(active[c1]);
	Link inv_lnk(c0, lnk.m.inverse());
	FOR(lf, 6) {
		Volume::Facet f0 = c0.facet(lf);
		Volume::Facet f1 = direction_to_face(c1, lnk.m * face_to_direction(f0));
		if (c0 == link[f0].q_id) return false;
		if (c1 == link[f1].q_id) return false;
	}
	return true;
}


void ReadOnlyMeshExtract3d::cancel(Volume::Facet seed_f) {
	auto lnk = link[seed_f];
	um_assert(lnk.active());
	Volume::Cell c0 = seed_f.cell();
	Volume::Cell c1(hex, lnk.q_id);
	um_assert(c0 != c1);
	um_assert(active[c0]);
	um_assert(active[c1]);
	Link inv_lnk(c0, lnk.m.inverse());
	FOR(lf, 6) {
		Volume::Facet f0 = c0.facet(lf);
		Volume::Facet f1 = direction_to_face(c1, lnk.m * face_to_direction(f0));
		if (c1 == link[f0].q_id) continue;
		if (c0 == link[f0].q_id) continue;

		Volume::Facet opp0 = link_opposite_facet(f0);
		Volume::Facet opp1 = link_opposite_facet(f1);

		if (opp1.active()) link[opp1] = link[opp1].followed_by(inv_lnk).followed_by(link[f0]);
		if (opp0.active()) link[opp0] = link[opp0].followed_by(lnk).followed_by(link[f1]);

		//link[f0] = lnk;
		//link[f1] = inv_lnk;

	}
	active[c0] = false;
	active[c1] = false;
}



void ReadOnlyMeshExtract3d::untangle() {
	//check that links are symmetric
	for (auto f : hex.iter_facets()) if (link[f].active()) um_assert(link_opposite_facet(link_opposite_facet(f)) == f);
	for (auto f : hex.iter_facets()) if (link[f].active()) um_assert(link_opposite_facet(f) != f);

	for (bool loop_only : {true, false}) {
		bool nothing_done;
		do {
			nothing_done = true;
			for (auto f : hex.iter_facets()) {
				if (!active[f.cell()]) continue;
				um_assert(f.active());
				if (!link[f].active()) continue;
				if (!active[link[f].q_id]) continue;
				auto opp = link_opposite_facet(f);
				if (loop_only && (emb[f.halfedge(0).opposite_f().facet()].f_id != -1 || emb[opp.halfedge(0).opposite_f().facet()].f_id != -1)) continue;

				if (reverted[f.cell()] == reverted[link[f].q_id]) continue;
				if (can_cancel(f))
				{
					cancel(f);
					nothing_done = false;
				}
			}
		} while (!nothing_done);
	}
	// show lasting revert
	CellAttribute<bool> badcell(hex);
	for (auto c : hex.iter_cells()) badcell[c] = active[c] && reverted[c];
	for (auto c : hex.iter_cells()) if (badcell[c]) Trace::alert("Reverted elements detected after untangling");
	Drop(hex, badcell)._skip_value(false).apply("BAD CELL");


	for (auto f : hex.iter_facets()) if (link[f].active() && active[f.cell()]) um_assert(link_opposite_facet(link_opposite_facet(f)) == f);
	for (auto f : hex.iter_facets()) if (link[f].active()) um_assert(link_opposite_facet(f) != f);
	for (auto f : hex.iter_facets()) if (link[f].active() && active[f.cell()]) um_assert(active[link_opposite_facet(f).cell()]);

}

// 
//
// void ReadOnlyMeshExtract3d::place_vertices() {
//
//	for (auto h : quad.iter_halfedges()) {
//
//		Surface::Facet f(tri, f_tri[h.facet()]);
//		vec2 pos = quad_center_U[h.facet()];
//
//		TransitionFunction tf;
//		Surface::Halfedge input_halfedge(tri, -1);
//
//		while (true) {
//			Triangle2 tr_U(uvw[f.halfedge(0)], uvw[f.halfedge(1)], uvw[f.halfedge(2)]);
//			FOR(lv, 3) tr_U[lv] = tf.apply(tr_U[lv]);
//
//			double f_orient = (tr_U.signed_area() > 0) ? 1 : (-1);
//
//			if (f_orient <= 0) break;
//
//			vec2 running_dir = .5 * (AxisDirection(h).vector() + AxisDirection(h.next()).vector());
//			vec2 obj = quad_center_U[h.facet()] + running_dir;
//
//
//			{
//				vec3 bc = tr_U.bary_coords(obj);
//				if (bc[0] > 0 && bc[1] > 0 && bc[2] > 0) 
//				{
//					h.to().pos() = bc[0] * f.vertex(0).pos() + bc[1] * f.vertex(1).pos() + bc[2] * f.vertex(2).pos();
//					break;
//				}
//			}
//
//			Surface::Facet next_f(tri, -1);
//			for (auto cir : f.iter_halfedges()) {
//				if (cir == input_halfedge) continue;
//				vec2 A = uvw[cir];
//				vec2 B = uvw[cir.next()];
//				A = tf.apply(A);
//				B = tf.apply(B);
//				vec2 I(0, 0);
//				if ((A - obj).norm2() < 1e-10) I = A;
//				if ((B - obj).norm2() < 1e-10) I = B;
//				if (I.norm2() > 0) {
//					vec3 bc = tr_U.bary_coords(I);
//					h.to().pos() = bc[0] * f.vertex(0).pos() + bc[1] * f.vertex(1).pos() + bc[2] * f.vertex(2).pos();
//					break;
//				}
//
//				double wB = Triangle2(pos, obj, A).signed_area();
//				double wA = -Triangle2(pos, obj, B).signed_area();
//				if (wA * wB <= 1e-20) continue;
//				I = (wB * B + wA * A) / (wB + wA);
//
//
//				if (!input_halfedge.active() && running_dir * (I - quad_center_U[h.facet()]) < 0) continue;
//				pos = I;
//
//
//				if (cir.opposite().active()) {
//
//					input_halfedge = cir.opposite();
//					next_f = cir.opposite().facet();
//					tf = tf.apply(TransitionFunction(cir, uvw));
//				}
//				h.to().pos() = (wB * cir.to().pos() + wA * cir.from().pos()) / (wB + wA);
//				break;
//			}
//			if (next_f.active()) f = next_f;
//			else break;
//		}
//	}
//}




// get projected dimension, position and anisotropy
std::tuple<int, vec3, mat3x3> get_fitting_constraint(std::vector<vec3>& pt, std::vector<vec3>& n) {
	mat3x3 AtA;
	vec3 AtB;

	vec3 G;
	FOR(p, pt.size()) G += pt[p];
	G /= double(pt.size());

	auto add_plane = [&](vec3 P, vec3 N) {
		FOR(i, 3)       FOR(j, 3)AtA[i][j] += N[i] * N[j];
		FOR(i, 3)       AtB[i] += P * N * N[i];
		};

	FOR(p, pt.size())   add_plane(pt[p], n[p]);
	FOR(i, 3)           AtA[i][i] += .00001;
	FOR(i, 3)           AtB[i] += .00001 * G[i];

	auto [eigen_val, eigen_vect] = eigendecompose_symmetric(AtA);

	std::tie(eigen_val, eigen_vect) = eigendecompose_symmetric(AtA);
	int dim = 0;
	if (eigen_val[1] < .5 * eigen_val[0]) {
		dim = 2;
		vec3 normal = eigen_vect.col(0);// or column ? idk
		vec3 x = cross(vec3(1, 0, 0), normal);
		if (x.norm2() < .1) x = cross(vec3(0, 1, 0), normal);
		x.normalize();
		vec3 y = cross(normal, x);
		add_plane(G, 10 * x);
		add_plane(G, 10 * y);
		std::tie(eigen_val, eigen_vect) = eigendecompose_symmetric(AtA);
	}
	else
		if (eigen_val[2] < .5 * eigen_val[0]) {
			dim = 1;
			vec3 x = cross(eigen_vect.col(0), eigen_vect.col(1));
			add_plane(G, 10 * x);
			std::tie(eigen_val, eigen_vect) = eigendecompose_symmetric(AtA);
		}
	mat3x3 inv_val;
	FOR(i, 3) inv_val[i][i] = 1. / eigen_val[i];
	mat3x3 eigen_inv = eigen_vect * inv_val * eigen_vect.transpose();

	return { dim,eigen_inv * AtB,eigen_vect };

}




void ReadOnlyMeshExtract3d::convert_to_hexes(Hexahedra& hexout, CellFacetAttribute<int>& emb_out) {

	// remove inactive hexes... an update links accordingly
	{
		std::vector<int> old2new(hex.ncells());
		CellAttribute<int> new2old(hex);
		for (auto c : hex.iter_cells()) new2old[c] = c;
		std::vector<bool> to_kill(hex.ncells());
		for (auto c : hex.iter_cells()) to_kill[c] = !active[c];
		hex.delete_cells(to_kill);
		hex.delete_isolated_vertices();


		for (auto c : hex.iter_cells()) old2new[new2old[c]] = c;
		for (auto f : hex.iter_facets()) if (link[f].q_id >= 0) link[f].q_id = old2new[link[f].q_id];

	}

	// Merge vertices
	DisjointSet ds(hex.nverts());
	for (auto h : hex.iter_halfedges()) {
		auto opp = link_opposite(h);
		if (!opp.active()) continue;
		ds.merge(h.to(), opp.from());
		ds.merge(h.from(), opp.to());
	}

	//PointAttribute<int> id2grp(hex);
	//int nverts = ds.get_sets_id(id2grp.ptr->data);
	//Drop(hex.points, id2grp).apply("grps");
	std::vector<int> id2grp;
	int nverts = ds.get_sets_id(id2grp);
	hexout.create_cells(hex.ncells());
	hexout.points.create_points(nverts);

	PointAttribute<int> grp_size(hexout, 0);
	FOR(v, hex.nverts()) {
		grp_size[id2grp[v]]++;
		hexout.points[id2grp[v]] += hex.points[v];
	}
	FOR(v, hexout.nverts())hexout.points[v] /= grp_size[v];
	FOR(c, hex.ncells()) FOR(lv, 8) hexout.vert(c, lv) = id2grp[8 * c + lv];
	hexout.connect();
	for (auto f : hexout.iter_facets()) emb_out[f] = emb[f].f_id;

	DropVolume(hexout).apply("hex_no_boundary_fit");
	ToolBox(hexout).drop_boundary();

	for (auto v : hex.iter_vertices()) v.pos() = .8 * v.pos() + .2 * hexout.points[id2grp[v]];
	DropVolume(hex).apply("hex"); return;



	// project boundary vertices on boundary
	std::map<int, std::vector<vec3> > boundaries;
	for (auto f : hexout.iter_facets()) {
		if (!active[f.cell()])continue;
		if (link[f].active()) continue;
		if (!f.on_boundary()) continue;
		if (emb[f].f_id == -1)  continue;
		for (auto h : f.iter_halfedges()) {
			Volume::Vertex v = h.from();
			auto it = boundaries.find(v);
			if (it == boundaries.end()) boundaries[v] = std::vector<vec3>();
			boundaries[v].push_back(emb[f].pos);
			boundaries[v].push_back(Triangle3(Volume::Facet(tet, emb[f].f_id)).normal());
		}
	}

	{// DEBUG OUTPUT START
		PointAttribute<int> nb(hexout, -1);
		for (auto v : hexout.iter_vertices()) {
			auto it = boundaries.find(v);
			if (it != boundaries.end()) nb[v] = it->second.size();
		}
		Drop(hexout.points, nb).apply("nb_adj_faces");

		CellFacetAttribute<int> embf(hexout, -1);
		for (auto f : hexout.iter_facets()) if (active[f.cell()] && !link[f].active() && f.on_boundary())
			embf[f] = emb[f].f_id;
		Drop(hexout, embf).apply("embf");
	}// DEBUG OUTPUT END

	for (auto v : hexout.iter_vertices()) {
		auto it = boundaries.find(v);
		if (it == boundaries.end()) continue;
		if (it->second.size() < 6) {
			Trace::alert("Boundary vertex with less than 3 boundary facets");
			continue;
		}
		std::vector<vec3> pt, n;
		for (auto data : it->second) if (pt.size() == n.size()) pt.push_back(data); else n.push_back(data);
		auto [dim, pos, ani] = get_fitting_constraint(pt, n);
		v.pos() = pos;
	}


}




void ReadOnlyMeshExtract3d::apply(Hexahedra& hex_out, CellFacetAttribute<int>& emb_out) {




	Trace::step("Pre-process");		pre_process();
	//return;
	Trace::step("Generate hexes");		generate_hexes();
	//show_hex_network(true, false);
	Trace::step("Connect hexes");		connect_hexes();
	//smooth_hex_network();
	//show_hex_network();
	Trace::step("Propagate boundary");	propagate_boundary();
	//sanity_check();

	show_hex_network(true, true);
	Trace::step("Untangle");			untangle();
	show_hex_network(true, true);
	//sanity_check();

	//Trace::step("Place vertices");	place_vertices();
	//sanity_check();	 

	Trace::step("Export hexes");
	convert_to_hexes(hex_out, emb_out);
}