#ifndef SDF__H__
#define	SDF__H__

#include <ultimaille/all.h>
#include <framework/trace.h>


#include "toolbox.h"
#include "drop_attribute.h"

#include "mesh_geom.h"

struct SDF {
	virtual ~SDF() {  }
	virtual double dist(vec3 P) = 0;
	virtual vec3 proj(vec3 P) = 0;
	virtual vec3 grad(vec3 P) = 0;
	void show(BBox3 box, int nb = 2000) {
		PolyLine pl;
		PointAttribute<double> dist_attr(pl.points, 0);
		box.dilate(.2 * (box.max - box.min).norm());
		pl.create_edges(nb);
		pl.points.create_points(2 * nb);
		FOR(sample, nb) {
			vec3 R;
			FOR(d, 3) R[d] = rand_range(box.min[d], box.max[d]);
			vec3 P = proj(R);
			FOR(i, 2) pl.vert(sample, i) = 2 * sample + i;
			pl.points[pl.vert(sample, 0)] = R;
			pl.points[pl.vert(sample, 1)] = P;
			dist_attr[pl.vert(sample, 0)] = dist(R);
		}

		DropPolyLine(pl).add(dist_attr, "dist")._active_point_attribute("dist").apply();
	}

};


/*
* hard coded explicit SDFs for testing
*/
struct CurvedCutSDF : public SDF {
	CurvedCutSDF(BBox3 box, int num = 0) : box(box), num(num) {}
	~CurvedCutSDF() {}
	double dist(vec3 P) {
		if (num == 5) {
			P = (P - .5 * (box.min + box.max)) / (box.max - box.min).norm();
			P[1] = 0;
			return P.norm() - .15;
		}

		if (num == 4) {
			P = (P - .5 * (box.min + box.max)) / (box.max - box.min).norm();
			return P.norm() - .15;
		}

		P = (P - box.min) / (box.max - box.min).norm();
		if (num == 0) return P.norm() - .5;
		if (num == 1) return (vec3(1, 0, 0) - P).norm() - .8;
		if (num == 2)  return P[0] - .3;
		if (num == 3)  return P[1] - .2;
		um_assert(false); return 0;
	};
	vec3 proj(vec3 P) {
		if (num == 5) {
			P = (P - .5 * (box.min + box.max)) / (box.max - box.min).norm();
			P = vec3(0, P[1], 0) + .15 * (P - vec3(0, P[1], 0)).normalized();
			P = P * (box.max - box.min).norm() + .5 * (box.min + box.max);
			return P;
		}

		if (num == 4) {
			P = (P - .5 * (box.min + box.max)) / (box.max - box.min).norm();
			P = .15 * P.normalized();
			P = P * (box.max - box.min).norm() + .5 * (box.min + box.max);
			return P;
		}

		P = (P - box.min) / (box.max - box.min).norm();
		if (num == 0) P = .5 * P.normalized();
		if (num == 1) P = vec3(1, 0, 0) - (.8 * (vec3(1, 0, 0) - P).normalized());
		if (num == 2) P[0] = .3;
		if (num == 3) P[1] = .2;
		P = P * (box.max - box.min).norm() + box.min;
		return P;
	}
	vec3 grad(vec3 P) {
		if (num == 5) {
			P = (P - .5 * (box.min + box.max)) / (box.max - box.min).norm();
			return (P - vec3(0, P[1], 0)).normalized();
		}

		if (num == 4) {
			P = (P - .5 * (box.min + box.max)) / (box.max - box.min).norm();
			return  P.normalized();
		}

		P = (P - box.min) / (box.max - box.min).norm();
		if (num == 0) return P.normalized();
		if (num == 1) return  (.8 * (vec3(1, 0, 0) - P).normalized());
		if (num == 2) return vec3(1,0,0);
		if (num == 3) return vec3(0, 1, 0);
		um_assert(false); return vec3(0,0,0);
	}

	int num;
	BBox3 box;
};


struct TrianglesSDF : public SDF {
	TrianglesSDF(Triangles& in) {
		ToolBox(m).copy_from(in);
		um_assert(m.nfacets() > 0);
		m.connect();
		npos = new NearestPointOnTriangles(m);
	}
	TrianglesSDF(Hexahedra& hex, CellFacetAttribute<int>& chart, int id) {
		ToolBox(m.points).copy_from(hex.points);
		FOR(f, hex.nfacets()) if (chart[f] == id) {
			int off_tri = m.create_facets(2);
			int quad2tri[2][3] = { {0,1,2},{0,2,3} };
			FOR(tri, 2)FOR(lv, 3)m.vert(off_tri + tri, lv) = hex.facet_vert(f / 6, f % 6, quad2tri[tri][lv]);
		}
		um_assert(m.nfacets() > 0);
		m.connect();
		npos = new NearestPointOnTriangles(m);
	}
	~TrianglesSDF() { delete npos; }
	double dist(vec3 P) {
		auto H = npos->request(P);
		return (P - H.pos) * Triangle3(H.f).normal();
	};
	vec3 grad(vec3 P) {
		auto H = npos->request(P);
		return  Triangle3(H.f).normal();
	};
	vec3 proj(vec3 P) {
		PointOnTriangles H = npos->request(P);
		return H;
	}
	Triangles m;
	NearestPointOnTriangles* npos;
};
#endif