#pragma once

#include <ultimaille/all.h>

using namespace UM;


struct HexBoundary {

	inline HexBoundary(Hexahedra& hex/*,bool duplicate_vertices=false*/) :hex(hex), quad_facet_(hex), quad(), hex_facet_(quad) {
 
		quad.points = hex.points;
		
		for (auto f_hex : hex.iter_facets()) if (!f_hex.opposite().active()) {
			int f_quad= quad.create_facets(1);
			hex_facet_[f_quad] = f_hex;
			quad_facet_[f_hex] = f_quad;
			for(int lv = 0; lv < 4; lv++) quad.vert(f_quad, lv) = f_hex.vertex(lv);
		}

		quad.connect();
	}

	inline HexBoundary(Hexahedra &hex, CellAttribute<bool> filter) : hex(hex), quad_facet_(hex), hex_facet_(quad), quad() {

		PointAttribute<bool> visible_point(hex, false);
		int n_verts = 0;
		for (auto c : hex.iter_corners()) {
			if (!filter[c.cell()] || visible_point[c.vertex()])
				continue;
			
			visible_point[c.vertex()] = true;
			++n_verts;
		}

		quad.points.create_points(n_verts);
		std::vector<int> hexverts_2_boundverts(hex.nverts(), -1);

		int i = 0;
		for (auto v : hex.iter_vertices()) {
			if (!visible_point[v])
				continue;

			quad.points[i] = hex.points[v];
			hexverts_2_boundverts[v] = i;
			++i;
		}

		for (auto f : hex.iter_facets()) {
			if (!filter[f.cell()])
				continue;

			int f_quad = quad.create_facets(1);
			hex_facet_[f_quad] = f;
			quad_facet_[f] = f_quad;

			for (int lv = 0; lv < 4; lv++) {
				quad.vert(f_quad, lv) = hexverts_2_boundverts[f.vertex(lv)];
			}
		}

		quad.connect();
	}

	inline void clear_surface() {
		quad.clear();
	}


	inline Volume::Facet		hex_facet(int f_quad) { return Volume::Facet(hex, hex_facet_[f_quad]); }
	inline Volume::Halfedge	hex_halfedge(int h_quad) { return Volume::Facet(hex, hex_facet_[h_quad/ 4]).halfedge(h_quad% 4); }

	inline Quads::Facet	quad_facet(int f_hex) { return Quads::Facet(quad, quad_facet_[f_hex]); }
	inline Quads::Halfedge	quad_halfedge(int h_hex) { return Quads::Facet(quad, quad_facet_[h_hex / 4]).halfedge(h_hex % 4); }


	inline void set_attribute_to_surface(CellFacetAttribute<int> &hex_attr, FacetAttribute<int> &quad_attr) {
		for (auto f : hex.iter_facets()) {
			if (!f.on_boundary())
				continue;
			
			quad_attr[quad_facet(f)] = hex_attr[f];
		}
	}

	inline void set_attribute_to_volume(FacetAttribute<int> &quad_attr, CellFacetAttribute<int> &hex_attr) {
		for (auto f : quad.iter_facets()) {			
			hex_attr[hex_facet(f)] = quad_attr[f];
		}
	}

	Hexahedra& hex;
	CellFacetAttribute<int> quad_facet_;

	Quads quad;
	FacetAttribute<int> hex_facet_;
};


struct TetBoundary {

	inline TetBoundary(Tetrahedra& tet) : tet(tet), tri_facet_(tet), tri(), tet_facet_(tri) {
		update();
	}

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

	inline void set_attribute_to_volume(FacetAttribute<int> &tri_attr, CellFacetAttribute<int> &tet_attr) {
		for (auto f : tri.iter_facets()) {			
			tet_attr[tet_facet(f)] = tri_attr[f];
		}
	}

	Tetrahedra& tet;
	CellFacetAttribute<int> tri_facet_;

	Triangles tri;
	FacetAttribute<int> tet_facet_;
};