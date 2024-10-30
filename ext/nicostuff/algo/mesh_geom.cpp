#include "mesh_geom.h"

#include "toolbox.h"
#include "drop_attribute.h"


using namespace UM;




void test_PointOnTriangles(Triangles& m) {
	PointOnTriangles P(Surface::Facet(m, 0), Triangle3(Surface::Facet(m, 0)).bary_verts());
	PolyLine pl;
	int nsteps = 50000;
	pl.create_edges(nsteps);
	pl.points.create_points(nsteps + 1);
	vec3 dir(1, 1, 0);
	double ave = ToolBox(m).ave_edge_size();
	FOR(t, nsteps) {
		pl.vert(t, 0) = t;
		pl.vert(t, 1) = t + 1;
		pl.points[t] = P;
		P.walk(dir, .5 * ave);
		pl.points[t + 1] = P;// needed only for last step
	}
	DropPolyLine(pl).apply("walk");
};


NearestPointOnTriangles::NearestPointOnTriangles(Triangles& m) : m(m) {
	bboxes.resize(m.nfacets());
	for (auto f : m.iter_facets())
		for (auto h : f.iter_halfedges())
			bboxes[f].add(h.from().pos());
	hbbox.init(bboxes);
}
PointOnTriangles NearestPointOnTriangles::request(vec3 v, double radius) {
		if (radius == -1) radius = (hbbox.tree[0].max - hbbox.tree[0].min).norm() / 100.;

	Surface::Facet  f_id(m, 0);
	vec3 P(0, 0, 0);

	double min_dist = 1e20;
	while (true) {
		std::vector<int> prim;
		BBox3 request; request.add(v); request.dilate(radius);
		hbbox.intersect(request, prim);
		for (int id : prim) {
			Surface::Facet f(m, id);
			Triangle3 t= Triangle3(f);
			vec3 bc = t.bary_coords(v);
			// test if on edge
			FOR(lh, 3) {
				vec3 nearest = Segment3(f.halfedge(lh)).nearest_point(v);
				double dist = (nearest - v).norm();
				if (min_dist > dist) { min_dist = dist; f_id = f; P=nearest;}
			};
			//// test if on triangle
			if (bc[0] > 0 && bc[1] > 0 && bc[2] > 0) {
				vec3 nearest = bc[0] * t[0] + bc[1] * t[1] + bc[2] * t[2];
				double dist = (nearest - v).norm();
				if (min_dist > dist) { min_dist = dist; f_id = f; P = nearest; }
			}
		}
		if (radius > min_dist)
			return { f_id, P };
		else radius *= 2.;
	}
}

InsideSurface::InsideSurface(Triangles& m) :m(m) {
	bboxes.resize(m.nfacets());
	for (auto f : m.iter_facets())
		for (auto h : f.iter_halfedges())
			bboxes[f].add(h.from().pos().xy());

	hbbox.init(bboxes);
};
bool InsideSurface::request(vec3 P) {
	int sum = 0;
	std::vector<int> prim;
	BBox2 request; request.add(P.xy()); request.dilate(.01);
	hbbox.intersect(request, prim);
	for (int id : prim) {
		Surface::Facet f(m, id);
		Triangle2 t2d = { f.vertex(0).pos().xy() ,f.vertex(1).pos().xy() ,f.vertex(2).pos().xy() };
		vec3 bc = t2d.bary_coords(P.xy());
		if (std::abs(t2d.signed_area()) < 1e-10) continue;
		if (bc[0] < 0 || bc[1] < 0 || bc[2] < 0) continue;

		Triangle3 t3d = { f.vertex(0).pos() ,f.vertex(1).pos() ,f.vertex(2).pos() };
		vec3 I = bc[0] * t3d[0] + bc[1] * t3d[1] + bc[2] * t3d[2];
		if (I[2] - P[2] < 0) continue;
		Tetrahedron tet( t3d[0], t3d[1], t3d[2], P );
		if (tet.volume() > 0) sum++; else sum--;
	}
	return sum != 0;
}









void test_NearestPointOnTriangles(Triangles& m) {
	NearestPointOnTriangles npos(m);
	PolyLine pl;
	auto box = ToolBox(m.points).bbox();
	box.dilate(.1 * (box.max - box.min).norm());

	int nb = 100000;
	pl.create_edges(nb);
	pl.points.create_points(2 * nb);
	FOR(sample, nb) {
		vec3 R;
		FOR(d, 3) R[d] = rand_range(box.min[d], box.max[d]);
		vec3 P = npos.request(R);

		FOR(i, 2) pl.vert(sample, i) = 2 * sample + i;
		pl.points[pl.vert(sample, 0)] = R;
		pl.points[pl.vert(sample, 1)] = P;
	}
	PointAttribute<int> start(pl.points);
	FOR(i, pl.nverts()) start[i] = i % 2;
	DropPolyLine(pl).add(start, "start")._active_point_attribute("start").apply();
}

void test_InsideSurface(Triangles& m) {
	InsideSurface is_in(m);
	PointSet pts;
	auto box = ToolBox(m.points).bbox();
	box.dilate(.1 * (box.max - box.min).norm());

	int nb = 100000;
	pts.create_points(nb);
	PointAttribute<int> in(pts);
	FOR(i, nb) {
		FOR(d, 3) pts[i][d] = rand_range(box.min[d], box.max[d]);
		in[i] = is_in.request(pts[i]);
	}
	Drop(pts, in).apply();

}
