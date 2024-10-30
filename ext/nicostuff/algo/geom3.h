#ifndef GEOM3__H__
#define  GEOM3__H__

#include <ultimaille/all.h>
#include <basic.h>

#include "geom2.h"

using namespace UM;




	namespace Intersect {

		bool triangle_line(Triangle3& tr, vec3 P, vec3 dir, vec3& bc);

		struct TetrahedraSegment {
			TetrahedraSegment(Tetrahedra& tet);
			bool inside(vec3 v);
			bool edge_intersects(vec3 org, vec3 dest);
			std::tuple<vec3, vec3> edge_intersection(vec3 org, vec3 dest);

			std::vector<BBox3> bboxes;
			HBoxes<3> hbbox;
			Tetrahedra& tet;
		};
	};

namespace Geom3{
	inline double vector_angle(vec3 v0, vec3 v1) { return atan2(cross(v0, v1).norm(), v0 * v1); }

	inline double dihedral_angle(Surface::Halfedge h){
		auto opp = h.opposite();
		um_assert(opp.active());
		vec3 n0,n1;
		if (h.facet().size()==3) n0= Triangle3(h.facet()).normal();
		else if (h.facet().size()==4) n0= Quad3(h.facet()).normal();
		else um_assert(false);

		if (h.facet().size()==3) n1= Triangle3(opp.facet()).normal();
		else if (h.facet().size()==4) n1= Quad3(opp.facet()).normal();
		else um_assert(false);
		return vector_angle(n0,n1);
	}



};

#endif