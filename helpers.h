
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

	std::vector<int> get_h_layers(Hexahedra &hex);

	std::vector<UM::vec3> get_layer_stack_visu(Hexahedra &hex, Volume::Halfedge &selected_he);
	std::vector<std::pair<int, int>> get_layer_stack(Hexahedra &hex, Volume::Halfedge &selected_he);
	void redefine_stack_layers(Hexahedra &hex, Volume::Halfedge &selected_he, int final_height);
	int get_facets_layers(UM::Hexahedra &hex, CellFacetAttribute<int> &layer);
	int get_layers(Hexahedra &hex, EdgeGraph &eg, EdgeAttribute<int> &layer);
	std::vector<int> get_cells_layer(Hexahedra &hex);
	std::vector<UM::vec3> remap_poly(std::vector<UM::vec3> poly, int n);

	void collapse(UM::Hexahedra &hex, std::vector<int> layer);

	void layer_from_halfedge(UM::Hexahedra &hex, int he, std::function<void(UM::Volume::Halfedge&)> f);
	void layer_along(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Facet&)> f);
	void loop_along(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge& /* cur_he */, bool/* on_border */)> f);

}