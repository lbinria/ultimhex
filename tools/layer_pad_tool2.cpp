#include "layer_pad_tool2.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"



bool LayerPad2::draw_object_properties() {
	ImGui::InputInt("Number of layer: ", &n_layers_requested);

	if(ImGui::Button("Redefine layers stack !")) {
		ctx.gui_mode = NewBlocPadding;
		return true;
	}

	return false;
}

void LayerPad2::draw_viewer_properties() {}

void LayerPad2::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
	gl_draw::draw_cells_overlay(ctx.mesh_, hovered_cells, colorMapInfo);
	gl_draw::draw_cell_facets_overlay(ctx.mesh_, selected_cell_facets, colorMapInfo, .5f);

}

void LayerPad2::hover_callback(double x, double y, int source) {
	
	if (!is_init_layers) {
		layers = helpers::get_h_layers(ctx.hex);
		is_init_layers = true;
	}

	if (ctx.is_cell_lfacet_hovered() && ctx.is_cell_edge_hovered()) {
		auto hovered_h = ctx.hovered_he(ctx.hex);
		hh = hovered_h;

		hovered_cells.clear();
		for (auto h : ctx.hex.iter_halfedges()) {
			if (layers[h] == layers[hovered_h])
				hovered_cells.push_back(Volume::Halfedge(ctx.hex, h).cell());
		}
	}
	
}

void LayerPad2::mouse_button_callback(int button, int action, int mods, int source) {

	if (!ctx.hex.connected() || !ctx.is_cell_hovered() || !ctx.is_cell_edge_hovered() || !ctx.is_cell_lfacet_hovered()) {
		return;
	}


	selected_h_idx = ctx.hovered_he(ctx.hex);
	Volume::Halfedge selected_he(ctx.hex, selected_h_idx);
	selected_cell_facets = helpers::get_layer_stack(ctx.hex, selected_he);

}

void LayerPad2::scroll_callback(double xoffset, double yoffset) {

}

void LayerPad2::validate_callback() {

	if (!ctx.is_cell_selected() || selected_h_idx < 0)
		return;


	std::cout << "n verts: " << ctx.hex.nverts() << std::endl;


	Volume::Halfedge selected_he(ctx.hex, selected_h_idx);
	helpers::redefine_stack_layers(ctx.hex, selected_he, n_layers_requested);
	um_bindings::geo_mesh_from_um_hex(ctx.hex, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	

	// Recompute layers
	ctx.hex.disconnect();
	ctx.hex.connect();
	// Reset
	clear();

	HexBoundary hb(ctx.hex);
	write_by_extension("bound_pour_voir_si_connectivite_pete.geogram", hb.quad);

	
	std::cout << "n verts: " << ctx.hex.nverts() << std::endl;
}

bool LayerPad2::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void LayerPad2::escape_callback() {
	// Clear tool
	clear();
}

