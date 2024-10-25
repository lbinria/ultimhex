#pragma once

#include <ultimaille/all.h>

using namespace UM;


struct HexBoundary {

	inline  HexBoundary(Hexahedra& hex/*,bool duplicate_vertices=false*/) :hex(hex), quad_facet_(hex), quad(), hex_facet_(quad) {
		// if (duplicate_vertices) 
		// 	ToolBox(quad.points).copy_from(hex.points);
		// else 
			quad.points = hex.points;
		
		for (auto f_hex : hex.iter_facets()) if (!f_hex.opposite().active()) {
			int f_quad= quad.create_facets(1);
			hex_facet_[f_quad] = f_hex;
			quad_facet_[f_hex] = f_quad;
			for(int lv = 0; lv < 4; lv++) quad.vert(f_quad, lv) = f_hex.vertex(lv);
		}
		quad.connect();
	}
	inline Volume::Facet		hex_facet(int f_quad) { return Volume::Facet(hex, hex_facet_[f_quad]); }
	inline Volume::Halfedge	hex_halfedge(int h_quad) { return Volume::Facet(hex, hex_facet_[h_quad/ 4]).halfedge(h_quad% 4); }

	inline Quads::Facet	quad_facet(int f_hex) { return Quads::Facet(quad, quad_facet_[f_hex]); }
	inline Quads::Halfedge	quad_halfedge(int h_hex) { return Quads::Facet(quad, quad_facet_[h_hex / 4]).halfedge(h_hex % 4); }

	Hexahedra& hex;
	CellFacetAttribute<int> quad_facet_;

	Quads quad;
	FacetAttribute<int> hex_facet_;
};

struct TetBoundary {

	inline TetBoundary(Tetrahedra& tet) : tet(tet), tri_facet_(tet), tri(), tet_facet_(tri) {
		// set(tet);
		update();
	}

	// inline void reset() {
	// 	for (auto f : tri.facets) {
	// 		tri_facet_[f] = -1;
	// 	}
	// 	for (int i = 0; i < tet.nfacets(); i++) {
	// 		tet_facet_[i] = -1;
	// 	}

	// 	tri.resize_attrs();

	// }

	inline void update() {
		
		tri.disconnect();
		tri.clear();

		// TODO not sure about this !
		tet.resize_attrs();
		tet_facet_.ptr.get()->resize(tet.nfacets());
		tri.resize_attrs();
		tri_facet_.ptr.get()->resize(tet.nfacets());

		tri.points = tet.points;
		for (auto f_tet : tet.iter_facets()) if (!f_tet.opposite().active()) {
			int f_tri = tri.create_facets(1);
			tet_facet_[f_tri] = f_tet;
			tri_facet_[f_tet] = f_tri;
			for (int lv = 0; lv < 3; lv++) tri.vert(f_tri, lv) = f_tet.vertex(lv);
		}

		tri.connect();
	}

	// inline void set(Tetrahedra &tet) {

	// 	tri.disconnect();

	// 	tri.points = tet.points;
	// 	for (auto f_tet : tet.iter_facets()) if (!f_tet.opposite().active()) {
	// 		int f_tri = tri.create_facets(1);
	// 		tet_facet_[f_tri] = f_tet;
	// 		tri_facet_[f_tet] = f_tri;
	// 		for (int lv = 0; lv < 3; lv++) tri.vert(f_tri, lv) = f_tet.vertex(lv);
	// 	}

	// 	tri.connect();
	// }

	inline Volume::Facet tet_facet(int f_tri) { return Volume::Facet(tet, tet_facet_[f_tri]); }

	inline Volume::Halfedge	tet_halfedge(int h_tri) { return Volume::Facet(tet, tet_facet_[h_tri / 3]).halfedge(h_tri % 3); }

	inline Triangles::Facet	tri_facet(int f_tet) { return Triangles::Facet(tri, tri_facet_[f_tet]); }
	
	inline Triangles::Halfedge	tri_halfedge(int h_tet) { return Triangles::Facet(tri, tri_facet_[h_tet / 3]).halfedge(h_tet % 3); }

	inline void set_attribute_to_surface(CellFacetAttribute<int> &tet_attr, FacetAttribute<int> &tri_attr) {
		for (auto f : tet.iter_facets()) {
			if (!f.on_boundary())
				continue;
			
			tri_attr[tri_facet(f)] = tet_attr[f];
		}
	}

	Tetrahedra& tet;
	CellFacetAttribute<int> tri_facet_;

	Triangles tri;
	FacetAttribute<int> tet_facet_;
};