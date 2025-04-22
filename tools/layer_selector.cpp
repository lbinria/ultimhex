#include "layer_selector.h"
#include "../geom_ultimaille_binding.h"
#include "../helpers.h"

void LayerSelector::init() {
	helpers::get_halfedge_layers(ctx.hex_bound->hex, layers);
}

void LayerSelector::hover() {

	Volume::Cell um_c(ctx.hex_bound->hex, ctx.hovered_cell);

	hovered_h = um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.hovered_edge, ctx.hovered_cell_lfacet));

	// Get all halfedges in layer of hovered halfedge
	hovered_halfedges.clear();
	for (auto h : ctx.hex_bound->hex.iter_halfedges()) {
		if (layers[h] == layers[hovered_h])
			hovered_halfedges.push_back(h);
	}
}

void LayerSelector::highlight_hovered_cells() {
	GEO::Attribute<int> cell_hovered_attr(
		ctx.mesh_.cells.attributes(), "cell_hovered"
	);

	// Remove last hovered cells
	for (auto hi : hovered_halfedges) {
		Volume::Halfedge h(ctx.hex_bound->hex, hi);

		if (cell_hovered_attr[h.cell()] == 1)
			cell_hovered_attr[h.cell()] = 0;
	}

	hover();

	// Set new hovered cells
	for (auto hi : hovered_halfedges) {
		Volume::Halfedge h(ctx.hex_bound->hex, hi);

		if (cell_hovered_attr[h.cell()] == 0)
			cell_hovered_attr[h.cell()] = 1;
	}
}

void LayerSelector::highlight_selected_cells() {
	GEO::Attribute<int> cell_hovered_attr(
		ctx.mesh_.cells.attributes(), "cell_hovered"
	);

	// Remove last selected
	for (auto hi : selected_halfedges) {
		Volume::Halfedge h(ctx.hex_bound->hex, hi);

		if (cell_hovered_attr[h.cell()] == 2)
			cell_hovered_attr[h.cell()] = 0;
	}

	select();

	for (auto hi : selected_halfedges) {
		Volume::Halfedge h(ctx.hex_bound->hex, hi);

		cell_hovered_attr[h.cell()] = 2;
	}
}

void LayerSelector::select() {
	selected_halfedges = hovered_halfedges;
	selected_layer = layers[hovered_h];
}

void LayerSelector::clear() {
	GEO::Attribute<int> cell_hovered_attr(
		ctx.mesh_.cells.attributes(), "cell_hovered"
	);
	
	// Remove last selected
	for (auto hi : selected_halfedges) {
		Volume::Halfedge h(ctx.hex_bound->hex, hi);

		if (cell_hovered_attr[h.cell()] == 2)
			cell_hovered_attr[h.cell()] = 0;
	}

	// Remove last hovered cells
	for (auto hi : hovered_halfedges) {
		Volume::Halfedge h(ctx.hex_bound->hex, hi);

		if (cell_hovered_attr[h.cell()] == 1)
			cell_hovered_attr[h.cell()] = 0;
	}


	hovered_halfedges.clear();
	selected_halfedges.clear();	
	hovered_h = -1;
	selected_layer = -1;
	layers.clear();
}