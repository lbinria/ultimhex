#include "path_constraint_padding_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"


void PathConstraintPaddingTool::slice() {



		




		// // Extract surface from selected facets
		// ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		// um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		// ctx.view.switch_to_surface_select_mode();

		// ctx.view.show_volume_ = true;
		// ctx.view.cells_shrink_ = 0.4f;
		// ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

}

void PathConstraintPaddingTool::reset() {
	clear();
	// Compute layers
	int nlayers = helpers::get_facets_layers(ctx.hex_bound->hex, layers);
}

bool PathConstraintPaddingTool::draw_object_properties() {

	if (ImGui::Button("Init##path_constraint_padding_tool_init")) {
		reset();
		// // Get layers
		// int nlayers = helpers::get_facets_layers(ctx.hex_bound->hex, layers);
		// selected_layers.resize(nlayers, -1);

		ctx.view.switch_to_volume_select_mode();
		
		ctx.view.cells_shrink_ = 0.5f;
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

		ctx.gui_mode = PathConstraintPadding;
		return true;
	}

	if (ctx.gui_mode != PathConstraintPadding && ctx.switch_mode != PathConstraintPadding)
		return false;

	ImGui::SliderFloat("Shrink", &ctx.view.cells_shrink_, 0, 1);

	return false;
}

void PathConstraintPaddingTool::draw_viewer_properties() {

}


void PathConstraintPaddingTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void PathConstraintPaddingTool::hover_callback(double x, double y, int source) {
	if (/*step > 0 ||*/ !ctx.is_facet_hovered())
		return;

}

void PathConstraintPaddingTool::mouse_button_callback(int button, int action, int mods, int source) {
	if (!ctx.is_cell_facet_hovered())
		return;

	// Get hovered layer
	int f = um_bindings::um_facet_index_from_geo_facet_index(ctx.hovered_cell_facet, 6);
	int l = layers[f];

	
	// if (ctx.left_mouse_pressed)
		selected_layers[l].push_back(f);
	// else 
		// selected_layers[l] = -1;

	// std::vector<int> layer_parts;

	{
	DisjointSet ds(ctx.hex_bound->hex.nfacets());

	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		// Skip facets not in a selected layer
		if (selected_layers[layers[f]].size() < 0)
			continue;

		// Check intersection with another selected layer
		for (auto h : f.iter_halfedges()) {
			bool intersect = false;
			for (auto eh : h.iter_CCW_around_edge()) {
				if (eh.facet() != f && selected_layers[layers[eh.facet()]].size() > 0 && layers[eh.facet()] != layers[f]) {
					intersect = true;
					break;
				}
			}

			if (!intersect && h.opposite_f().opposite_c().active())
				ds.merge(f, h.opposite_f().opposite_c().opposite_f().facet());
		}
		
	}

	int n_parts = ds.get_sets_id(layer_parts);
	std::cout << "n_parts: " << n_parts << std::endl;
	}

	

	// Only keep layer parts that contains the selected facets
	// std::vector<bool> selected_layer_parts(layer_parts.size(), false);
	CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		if (selected_layers[layers[f]].size() < 0)
			continue;

		auto layer_facets = selected_layers[layers[f]];

		for (auto lf : layer_facets) {
			if (layer_parts[f] == layer_parts[lf]) {
				selected[f] = true;
			}
		}
	}



	// CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
	// for (auto f : ctx.hex_bound->hex.iter_facets()) {

	// }


	Quads q;
	q.points.create_points(ctx.hex_bound->hex.nverts());
	q.create_facets(ctx.hex_bound->hex.nfacets());

	for (auto v : ctx.hex_bound->hex.iter_vertices()) {
		q.points[v] = ctx.hex_bound->hex.points[v];
	}

	// CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
	FacetAttribute<int> colors(q, -1);
	for (auto f: ctx.hex_bound->hex.iter_facets()) {
		if (selected_layers[layers[f]].size() < 0)
			continue;

		q.vert(f, 0) = f.vertex(0);
		q.vert(f, 1) = f.vertex(1);
		q.vert(f, 2) = f.vertex(2);
		q.vert(f, 3) = f.vertex(3);

		colors[f] = layer_parts[f];

		// selected[f] = true;
	}

	// colors.ptr->data = layer_parts;

	write_by_extension("padus.geogram", q, {{}, {{"f", colors.ptr}}, {}});


	// Update view
	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
	ctx.view.switch_to_volume_select_mode();
	ctx.view.show_surface_ = true;
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);


}

void PathConstraintPaddingTool::scroll_callback(double xoffset, double yoffset) {}

void PathConstraintPaddingTool::validate_callback() {

	CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		if (selected_layers[layers[f]].size() < 0)
			continue;

		auto layer_facets = selected_layers[layers[f]];

		for (auto lf : layer_facets) {
			if (layer_parts[f] == layer_parts[lf]) {
				selected[f] = true;
			}
		}
	}

	ctx.hex_bound->clear_surface();
	BenjaminAPI::pad(ctx.hex_bound->hex, selected, *ctx.emb_attr);

	// Reconstruct
	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);

	// Update view
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	ctx.view.switch_to_volume_select_mode();
	ctx.mesh_gfx_.unset_filters();
	ctx.view.show_vertices_ = false;

	reset();

}

void PathConstraintPaddingTool::escape_callback() {
	// Just clear selected layers
	selected_layers.clear();
	layer_parts.clear();
	
	// Reconstruct
	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);

	// Update view
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	ctx.view.switch_to_volume_select_mode();
	ctx.mesh_gfx_.unset_filters();
	ctx.view.show_vertices_ = false;
}

bool PathConstraintPaddingTool::is_compatible() { 
	return true;
}