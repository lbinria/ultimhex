#ifndef MESH_GEOM__H__
#define  MESH_GEOM__H__
#include <ultimaille/all.h>
#include <basic.h>

#include <toolbox_triangles.h>

using namespace UM;




struct PointOnTriangles {
	PointOnTriangles(Surface::Facet f, vec3 P) : f(f), pos(P) {}
	
	vec3&  project(){
		vec3 &P = pos;
		vec3 bc = Triangle3(f).bary_coords(P);
		FOR(lv,3) bc[lv] = std::max(.001,std::min(.999,bc[lv]));
		double sum = bc[0]+bc[1]+bc[2];
		FOR(lv,3) bc[lv] /= sum;
		P = bc[0]*f.vertex(0).pos()+bc[1]*f.vertex(1).pos()+bc[2]*f.vertex(2).pos();
		return P;
	}

	bool walk(vec3& dir, double step,CornerAttribute<bool>* feature_edge = NULL) {
		vec3& P = pos;
		vec3 n = Triangle3(f).normal();
		dir.normalize();
		// make sure that dir is normalized and in tangent space
		dir = dir - n * dir * n;				
		if (dir.norm2() < 1e-15) return false;
		step*=dir.norm();
		dir.normalize();

		while (step > 1e-10) {
			vec3 O = P + step * dir;
			// find the edge to escape the triangle
			double best_lambda = 0;
			Surface::Halfedge best_h = f.halfedge(0);
			best_h = -1;
			FOR(lh, 3) {
				auto h = f.halfedge(lh);
				vec3 edge = h.to().pos() - h.from().pos();
				double area[2] = {
					n * cross(edge,P - h.from().pos()),
					n * cross(O - h.from().pos(),edge),
				};
				if (area[1] + area[0]<1e-10 )					continue;	// dir is not going to cross the edge

				double lambda = area[0] / (area[0] + area[1]);
				if (best_h.active() && best_lambda < lambda)	continue;	// cross edge earlier than the other					
				best_h = h;
				best_lambda = lambda;
			}
			if (!best_h.active()) return false;
			// is it solved ?
			if (best_lambda > 1) { P = O; return true; }
			
			P = P + best_lambda * (O - P);
			auto opp = best_h.opposite();
			if (!opp.active()) return false; // reach a boundary
			if (feature_edge!=NULL && (*feature_edge)[best_h]) return false; // reach a boundary

			// parallel transort the dir vector to the opposite facet
			vec3 e_dir = (best_h.to().pos() - best_h.from().pos()).normalized();
			vec3 perp_f = cross(e_dir, n);
			vec3 n_opp = Triangle3(opp.facet()).normal();
			vec3 perp_opp = cross(e_dir, n_opp);
			n = n_opp;
			dir = dir * e_dir * e_dir + dir * perp_f * perp_opp;
			f = opp.facet();
			step = (1. - best_lambda) * step;
		}
		return false;
	}
	operator vec3()  { return pos; }

	Surface::Facet f;
	vec3 pos;
};

void test_PointOnTriangles(Triangles& m) ;








struct PointOnPolyLine :vec3 {
  PointOnPolyLine(PolyLine::Edge e, vec3 P) : e(e), vec3(P) {}

  vec3& project() {
    vec3& P = *this;
    double bc0 = Segment3(e).bary_coords(P);
    bc0 = std::max(0.001, std::min(.999, bc0));
    P = bc0 * e.from().pos() + (1. - bc0) * e.to().pos();
	return P;
  }


  bool walk(vec3& dir, double step) {
    //return true;
    vec3& P = *this;
    um_assert(e.opposite().active());//.warning: may have problems woth boundary
    if (Segment3(e).vector() * dir < 0) e = e.opposite();

    step *= Segment3(e).vector().normalized() * dir;

    // make sure that dir is normalized and in tangent space
    while (step > 1e-10) {
      um_assert(e.opposite().active());//.warning: may have problems woth boundary

      double dist_to_edge_extremity = (e.to().pos() - P).norm();
      if (dist_to_edge_extremity > step) {
        P += step * Segment3(e).vector().normalized();
        return true;
      }

      P = e.to().pos();
      step = step - dist_to_edge_extremity;
      auto next = e;
      for (auto cir : e.to().iter_edges()) if (cir != e.opposite()) {
        if (next == e) next = cir;
        else return false;
      }
      if (e == next)return false; // dandling feature 
      e = next;
    }
	return false;
  }

  PolyLine::Edge e;
};



struct PointInMeshData {
  int dimension = -1;
  int id = -1;
  vec3 P;
};




struct PointInMesh {
  PointInMesh(Triangles* tri, PolyLine* pl, PointInMeshData& data) : tri(tri), pl(pl), data(data) {}

  void project() {
    if (data.dimension == 1) {
      PointOnPolyLine pt(PolyLine::Edge(*pl, data.id), data.P);
      pt.project();
      data.P = pt;
    }
    if (data.dimension == 2) {
      PointOnTriangles pt(Surface::Facet(*tri, data.id), data.P);
      pt.project();
      data.P = pt;
    }
  }

  std::array<vec3, 3> constrained_directions() {
    if (data.dimension == 2) return { Triangle3(Surface::Facet(*tri,data.id)).normal(),vec3(0,0,0),vec3(0,0,0) };
    if (data.dimension == 1) {
      vec3 dir = Segment3(PolyLine::Edge(*pl, data.id)).vector().normalized();
      vec3 x(1, 0, 0);
      if (std::abs(dir * x) > .8) x = vec3(0, 1, 0);
      vec3 y = cross(dir, x).normalized();
      x = cross(dir, y).normalized();
      return { x,y,vec3(0,0,0) };
    }
    if (data.dimension == 0) return { vec3(1,0,0),vec3(0,1,0),vec3(0,0,1) };
    return { vec3(0,0,0),vec3(0,0,0),vec3(0,0,0) };
  }

  bool walk(vec3& dir, double step, CornerAttribute<bool>* feature_edge = NULL) {
    if (data.dimension == 0) return true;
    if (data.dimension == 1) {
      PointOnPolyLine pt(PolyLine::Edge(*pl, data.id), data.P);
      if (!pt.walk(dir, step)) return false;
      data.P = pt;
      data.id = pt.e;
      return true;
    }
    if (data.dimension == 2) {
      PointOnTriangles pt(Surface::Facet(*tri, data.id), data.P);
      if (!pt.walk(dir, step, feature_edge)) return false;
      data.P = pt;
      data.id = pt.f;
      return true;
    }
	return false;
  }
  vec3& pos() { return data.P; }
  PointInMeshData& data;
  Triangles* tri;
  PolyLine* pl;
};





struct NearestPointOnTriangles {
	Triangles& m;
	std::vector<BBox3> bboxes;
	HBoxes<3> hbbox;

	NearestPointOnTriangles(Triangles& m);
	PointOnTriangles request(vec3 v, double radius = -1);
};


struct TrianglesSegmentIntersections {
	TrianglesSegmentIntersections(Triangles& tri) : m(tri) {
		bboxes.resize(m.nfacets());
		for (auto f : m.iter_facets())
			for (auto h : f.iter_halfedges())
				bboxes[f].add(h.from().pos());
		hbbox.init(bboxes);
	}

	std::vector<PointOnTriangles> request(Segment3 seg) {
		std::vector<PointOnTriangles> res;
		std::vector<int> prim;
		BBox3 request;
		FOR(lv, 2) request.add(seg[lv]);
		hbbox.intersect(request, prim);
		for (int id : prim) {
			Surface::Facet f(m, id);
			vec3 bc;
			FOR(lh, 3) bc[lh] = Tetrahedron(f.vertex((lh + 1) % 3).pos(), f.vertex((lh + 2) % 3).pos(), seg[0], seg[1]).volume();
			double sum = 0;
			FOR(lh, 3) sum += bc[lh];
			FOR(lh, 3) bc[lh] /= sum;
			if (bc[0] < 0 || bc[1] < 0 || bc[2] < 0) continue;
			//FOR(lh, 3) plop(bc[lh]);
			vec3 I = bc[0] * f.vertex(0).pos() + bc[1] * f.vertex(1).pos() + bc[2] * f.vertex(2).pos();
			if ((seg[0] - I) * (seg[1] - I) > 0) continue;

			res.push_back(PointOnTriangles(f, I));
		}
		//std::sort(res.begin(), res.end(), [&](const PointOnTriangles& a, const PointOnTriangles& b) {
		//	return (seg[0] - a).norm2() > (seg[0] - b).norm2();
		//	});
		return res;
	}


	std::vector<BBox3> bboxes;
	HBoxes<3> hbbox;
	Triangles& m;
};


/*
* WARNINGS:
*    no robust predicates here
*    assumes that the surface is triangulated
*/
struct InsideSurface {
	InsideSurface(Triangles& m) ;
	bool request(vec3 P);

	std::vector<BBox2> bboxes;
	HBoxes<2> hbbox;
	Triangles& m;
};









void test_NearestPointOnTriangles(Triangles& m);
void test_InsideSurface(Triangles& m);

#endif

