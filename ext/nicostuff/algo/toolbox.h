#ifndef TOOLBOX__H__
#define TOOLBOX__H__

#define REAL double
#include <third_party/triangle.h>

#include <ultimaille/all.h>
#include <basic.h>
#include <framework/trace.h>
#include <third_party/tetgen/tetgen.h>


#include "geom2.h"
#include "geom3.h"
#include <map>

#undef min//#include <windows.h>


#include "toolbox_pointset.h"
#include "toolbox_polyline.h"
#include "toolbox_triangles.h"

#include "toolbox_tetrahedra.h"

#include "toolbox_hexahedra.h"
using namespace UM;

template<>
struct ToolBox<Surface::Halfedge> {
	ToolBox(Surface::Halfedge h) :h(h) {  }
	double corner_angle() {
		return Geom3::vector_angle(h.to().pos() - h.from().pos(), h.prev().from().pos() - h.from().pos());
	}

	double dihedral_angle() {
		um_assert(h.opposite().active());
		return std::acos(std::min(1., std::max(-1., Triangle3(h.facet()).normal() * Triangle3(h.opposite().facet()).normal())));
	}

	Surface::Halfedge next_on_border() {
		for (auto& it : h.to().iter_halfedges()) if (!it.opposite().active()) return it;
		um_assert(false);
	}
	Surface::Halfedge h;
};

template<>
struct ToolBox<Surface::Vertex> {
	ToolBox( Surface::Vertex v) :v(v) {  }
	int valence() {
		int ret = 0;
		for (Surface::Halfedge cir : v.iter_halfedges()) {
			um_assert(cir.active());
			ret++;
		}
		return ret;
	}
	Surface::Vertex v;
};




template<>
struct ToolBox<Surface> {
	ToolBox(Surface& m) :m(m) {  }
	double ave_edge_size() {
		double sum = 0;
		FOR(f, m.nfacets()) FOR(lv, 3) sum += (m.points[m.vert(f, lv)] - m.points[m.vert(f, (lv + 1) % 3)]).norm();
		return sum / double(3 * m.nfacets());
	}

	vec3 barycenter(int f) {
		vec3 res(0, 0, 0);
		FOR(lv, m.facet_size(f)) res += m.points[m.vert(f, lv)];
		return res / double(m.facet_size(f));
	}

	Surface& m;
};

template<>
struct ToolBox<Quads> : public ToolBox<Surface> {
	ToolBox(Quads& m) :ToolBox<Surface>(m) {  }
};

template<>
struct ToolBox<Volume::Halfedge> {
	ToolBox(Volume::Halfedge h) :h(h) {  }
	double dihedral_angle() {
		return std::acos(std::min(1., std::max(-1., Triangle3(h.facet()).normal() * Triangle3(h.opposite_f().facet()).normal())));
	}
	double corner_angle() {
		return Geom3::vector_angle(h.to().pos() - h.from().pos(), h.prev().from().pos() - h.from().pos());
	}
	Volume::Halfedge h;
};





#endif