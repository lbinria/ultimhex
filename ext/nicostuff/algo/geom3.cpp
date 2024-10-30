#include <geom3.h>
namespace Intersect {


	bool triangle_line(Triangle3& tr, vec3 P, vec3 dir, vec3& bc) {
		vec3 n = tr.normal();
		double proj_n = dir * n;
		double w = .5;
		if (std::abs(proj_n) < 1e-15) return false;
		w = (tr[0] - P) * n / proj_n;
		if (w < 0) return false;
		vec3 I_uvw = P + w * dir;
		bc = tr.bary_coords(I_uvw);
		FOR(lv, 3) if (bc[lv] < -1e-10) return false;
		return true;
	}

	TetrahedraSegment::TetrahedraSegment(Tetrahedra& tet) : tet(tet) {
		tet.connect();
		bboxes.resize(tet.ncells());
		FOR(c, tet.ncells()) FOR(lv, 4) bboxes[c].add(tet.points[tet.vert(c, lv)]);
		hbbox.init(bboxes);
	}
	bool TetrahedraSegment::inside(vec3 v) {
		std::vector<int> prim;
		BBox3 request; request.add(v); request.dilate(.001);
		hbbox.intersect(request, prim);
		for (int c : prim) {
			Tetrahedron t(tet.points[tet.vert(c, 0)], tet.points[tet.vert(c, 1)], tet.points[tet.vert(c, 2)], tet.points[tet.vert(c, 3)]);
			vec4 bc = t.bary_coords(v);
			bool is_in = true;
			FOR(lv, 4) is_in = is_in && (bc[lv] > -1e-20);
			if (is_in) return true;
		}
		return false;
	}


	bool TetrahedraSegment::edge_intersects(vec3 org, vec3 dest) { return inside(org) != inside(dest); }


	std::tuple<vec3, vec3> TetrahedraSegment::edge_intersection(vec3 org, vec3 dest) {
		um_assert(edge_intersects(org, dest));
		std::vector<int> prim;
		BBox3 request; request.add(org);
		request.add(dest);
		request.dilate(.001);
		hbbox.intersect(request, prim);
		for (int c : prim) FOR(lf, 4) {
			if (tet.conn->oppf[4 * c + lf] != -1) continue;
			Triangle3 tr({ tet.points[tet.facet_vert(c,lf, 0)], tet.points[tet.facet_vert(c,lf, 1)], tet.points[tet.facet_vert(c,lf, 2)] });
			vec3 n = tr.normal();
			double proj_n = (dest - org) * n;
			double w = .5;
			if (std::abs(proj_n) > 1e-15) w = (tr[0] - org) * n / proj_n;

			if (w < 0 || w>1) continue;
			vec3 I = org + w * (dest - org);
			vec3 bc = tr.bary_coords(I);

			bool is_in = true;
			FOR(lv, 3) is_in = is_in && (bc[lv] > -.1);
			if (is_in && std::abs(n * (I - tr[0])) < 1e-5) return { I,n };

		}
		vec3 P;
		P = (org + dest) / 2.;
		FOR(it, 5) {
			if (edge_intersects(P, dest))	org = P;
			else							dest = P;
			P = 0.5 * (org + dest);
		}
		return { P,vec3(0,0,0) };


	}

};
