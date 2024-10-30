#include <ultimaille/all.h>

#include <algo/drop_glyph.h>



void drop_point(vec3 P, std::string name) {
	PointSet pts;
	pts.create_points(1);
	pts[0] = P;
	DropPointSet(pts).apply(name);
}
void drop_arrow(vec3 from, vec3 to, std::string name) {
	PolyLine pl;
	pl.create_edges(1);
	pl.vert(0, 0) = 0;
	pl.vert(0, 1) = 1;
	pl.points.create_points(2);
	pl.points[0] = from;
	pl.points[1] = to;
	DropPolyLineGeometry(pl).apply_arrow(name);
}

void drop_triangle(Triangle3 t, std::string name) {
	Triangles m;
	m.points.create_points(3);
	m.create_facets(1);
	FOR(lv, 3) m.points[lv] = t[lv];
	FOR(lv, 3) m.vert(0, lv) = lv;
	DropSurface(m).apply(name);
}

