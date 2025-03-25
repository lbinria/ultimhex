#include "hex_collapse_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"


bool HexCollapseTool::draw_object_properties() {

	if (ImGui::Button("Select hex layer to collapse")) {
		ctx.gui_mode = HexCollapse;
		reset();
		return true;
	}

	ImGui::Checkbox("Auto smooth#chk_auto_smooth_collapse_tool", &auto_smooth);

	return false;
}

void HexCollapseTool::compute_layers() {
	// Compute hex layers
	DisjointSet ds(ctx.hex_bound->hex.ncells() * 24);

	for (auto h : ctx.hex_bound->hex.iter_halfedges()) {

		auto opp = h.opposite_f().opposite_c();
		if (opp.active())
			ds.merge(h, opp.opposite_f().next().next());
			
		opp = h.opposite_c();
		if (opp.active()) {
			ds.merge(h, opp.opposite_f().next().next().opposite_f());
		}
		
	}

	ds.get_sets_id(layers);
}



void HexCollapseTool::draw_viewer_properties() {}

void HexCollapseTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
}

void HexCollapseTool::hover_callback(double x, double y, int source) {

	if (!ctx.is_cell_edge_hovered() || !ctx.is_cell_lfacet_hovered())
		return;


	Volume::Cell um_c(ctx.hex_bound->hex, ctx.hovered_cell);

	hovered_h = um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.hovered_edge, ctx.hovered_cell_lfacet));


	hovered_cells.clear();
	for (auto h : ctx.hex_bound->hex.iter_halfedges()) {
		if (layers[h] == layers[hovered_h])
			hovered_cells.push_back(h.cell());
	}

	GEO::Attribute<int> cell_hovered_attr(
		ctx.mesh_.cells.attributes(), "cell_hovered"
	);


	// Remove last hovered cells
	for (auto c : last_hovered_cells) {
		if (cell_hovered_attr[c] == 1)
			cell_hovered_attr[c] = 0;
	}

	last_hovered_cells = hovered_cells;

	// Set new hovered cells
	for (auto c : hovered_cells) {
		if (cell_hovered_attr[c] == 0)
			cell_hovered_attr[c] = 1;
	}

	ctx.view.change_mode(ViewBinding::Mode::Volume);
	ctx.view.attribute_subelements_ = GEO::MeshElementsFlags::MESH_CELLS;
	ctx.view.attribute_ = "cells.cell_hovered";
	ctx.view.attribute_name_ = "cell_hovered";
	ctx.view.attribute_min_ = 0;
	ctx.view.attribute_max_ = 2;
}

void HexCollapseTool::mouse_button_callback(int button, int action, int mods, int source) {

	GEO::Attribute<int> cell_hovered_attr(
		ctx.mesh_.cells.attributes(), "cell_hovered"
	);

	// Remove last selected
	for (auto c : selected_cells) {
		if (cell_hovered_attr[c] == 2)
			cell_hovered_attr[c] = 0;
	}

	selected_cells = hovered_cells;


	for (auto c : selected_cells) {
		cell_hovered_attr[c] = 2;
	}

	selected_layer = layers[hovered_h];
}

void HexCollapseTool::scroll_callback(double xoffset, double yoffset) {

}

void HexCollapseTool::validate_callback() {
	
	std::vector<int> selected_halfedges;

	for (auto h : ctx.hex_bound->hex.iter_halfedges()) {
		if (layers[h] != selected_layer)
			continue;

		selected_halfedges.push_back(h);

		auto nh = h.opposite_f().next().next();
		selected_halfedges.push_back(nh);
		nh = nh.opposite_f().next().next();
		selected_halfedges.push_back(nh);
		nh = nh.opposite_f().next().next();
		selected_halfedges.push_back(nh);

	}

	{
		// Visualize embedding before collapse
		FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
		ctx.hex_bound->set_attribute_to_surface(*ctx.emb_attr, surf_emb_attr);
		write_by_extension("before_collapse_emb.geogram", ctx.hex_bound->quad, {{}, {{"emb", surf_emb_attr.ptr}}, {}});
	}

	// Clear temporary because we'll delete points shared between hex / quad of hex bound
	ctx.hex_bound->clear_surface();
	helpers::collapse(ctx.hex, selected_halfedges, *ctx.emb_attr);

	if (auto_smooth)
		BenjaminAPI::smooth(ctx.hex_bound->hex, *ctx.emb_attr, ctx.tet_bound->tri, *ctx.tri_chart);

	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

	{
		// Visualize feedback
		// CellFacetAttribute<int> emb(ctx.hex_bound->hex, -1);
		// emb.ptr->data = ctx.emb;
		FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
		ctx.hex_bound->set_attribute_to_surface(*ctx.emb_attr, surf_emb_attr);
		write_by_extension("collapse_emb.geogram", ctx.hex_bound->quad, {{}, {{"emb", surf_emb_attr.ptr}}, {}});
	}
	reset();
}

bool HexCollapseTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void HexCollapseTool::escape_callback() {

}