#include "facet_selector.h"
#include "../geom_ultimaille_binding.h"

void FacetSelector::paint(unsigned int hovered_facet) {
	std::vector<int> hovered_facets;

	// Chart selection mode
	// if (select_mode == 1) {
	// 	for (auto f : ctx.hex_bound->quad.iter_facets()) {
	// 		if (patches[ctx.hex_bound->hex_facet(f)] == patches[ctx.hex_bound->hex_facet(ctx.hovered_facet)])
	// 			hovered_facets.push_back(f);
	// 	}
	// } else {
		hovered_facets.push_back(hovered_facet);
	// }

	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	// Reset last hovered facets as not hovered, if not selected
	for (auto last_hovered_f : last_hovered_facets) {
		if (last_hovered_f >= 0 && hovered_attr[last_hovered_f] == 1) {
			hovered_attr[last_hovered_f] = 0;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
		}
	}

	// Last hovered facet become current hovered facet
	last_hovered_facets = hovered_facets;

	for (auto f : hovered_facets) {
		// If not selected, facet is hovered
		if (hovered_attr[f] < 2) {
			hovered_attr[f] = 1;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 1;
		}

		if (ctx.left_mouse_pressed) {
			hovered_attr[f] = 2;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 2;
		} else if (ctx.right_mouse_pressed) {
			hovered_attr[f] = 0;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 0;
		}
	}

}

void FacetSelector::clear() {
	last_hovered_facets.clear();
}