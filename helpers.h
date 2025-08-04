
#pragma once

#include <ultimaille/all.h>
#include <set>
#include <map>

using namespace UM;

namespace helpers {

	struct EdgeCompare {
		bool operator() (const PolyLine::Edge& lhs, const PolyLine::Edge& rhs) const {
			int e1 = lhs; int e2 = rhs;
			return e1 < e2;
		}
	};
	struct FacetComparator {
		bool operator() (const Volume::Facet& lhs, const Volume::Facet& rhs) const {
			int e1 = lhs; int e2 = rhs;
			return e1 < e2;
		}
	};

	void puff(UM::Hexahedra &hex, UM::Volume::Halfedge start_h, std::vector<int> layer, CellFacetAttribute<bool> &selected);
	void cross_puff(UM::Hexahedra &hex, std::vector<int> layer, CellFacetAttribute<bool> &selected);
	
	std::vector<int> get_h_layers(Hexahedra &hex);

	std::vector<std::pair<int, int>> get_layer_stack_facets(Hexahedra &hex, Volume::Halfedge &selected_he);
	void redefine_stack_layers(Hexahedra &hex, Volume::Halfedge &selected_he, int final_height);
	int get_facets_layers(UM::Hexahedra &hex, CellFacetAttribute<int> &layer);
	int get_facets_layers(UM::Hexahedra &hex, std::vector<int> &layer);
	int get_layers(Hexahedra &hex, EdgeGraph &eg, EdgeAttribute<int> &layer);
	int get_halfedge_layers(Hexahedra &hex, std::vector<int> &layer);

	std::vector<int> get_cells_layer(Hexahedra &hex);
	std::vector<UM::vec3> remap_poly(std::vector<UM::vec3> poly, int n);

	void collapse(UM::Hexahedra &hex, std::vector<int> layer, CellFacetAttribute<int> &emb_attr);

	void layer_from_halfedge(UM::Hexahedra &hex, int he, std::function<void(UM::Volume::Halfedge&)> f);
	void layer_along(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Facet&)> f);
	void loop_along(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge& /* cur_he */, bool/* on_border */)> f);


	// Given a reference cell, a halfedge, return combinatorial coordinates of the 8 corners of the hex cell
	inline std::array<std::pair<int, int>, 8> get_hex_corners(Hexahedra &m, int hi) {
		Volume::Halfedge h(m, hi);
		auto c = h.cell();

		std::vector<int> facets = {
			h.prev().opposite_f().facet(), 
			h.next().opposite_f().facet()
		};

		std::array<std::pair<int, int>, 8> corners;
		for (int lc = 0; lc < 8; ++lc) {
			for (int i = 0; i < facets.size(); ++i) {
				for (int lv = 0; lv < 4; ++lv) {
					// auto f = facets[i];
					auto f = Volume::Facet(m, facets[i]);
					if (c.corner(lc) == f.corner(lv)) {
						corners[lc] = {i, lv};
						break;
					}
				}
			}
		}

		return corners;
	}

}