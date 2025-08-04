#include "hex_split_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"

#include <algorithm>


void HexSplitTool::reset() {
	clear();
	layer_selector.init();
}

bool HexSplitTool::draw_object_properties() {

	if (ImGui::Button("Select hex layer to split##btn_hex_split_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		
		ctx.gui_mode = HexSplit;
		reset();
		
		ctx.view.switch_to_volume_select_mode();
		ctx.view.cells_shrink_ = 0.0f;

		return true;
	}

	ImGui::Checkbox("Auto smooth##chk_auto_smooth_hex_split_tool", &auto_smooth);

	return false;
}

void HexSplitTool::draw_viewer_properties() {}

void HexSplitTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
}

void HexSplitTool::hover_callback(double x, double y, int source) {

	if (!ctx.is_cell_edge_hovered() || !ctx.is_cell_lfacet_hovered())
		return;

	layer_selector.highlight_hovered_cells();


}

void HexSplitTool::mouse_button_callback(int button, int action, int mods, int source) {

	layer_selector.highlight_selected_cells();

}

void HexSplitTool::scroll_callback(double xoffset, double yoffset) {

}

void HexSplitTool::validate_callback() {
	
	auto selected_halfedges = layer_selector.get_selected_halfedges();

	HalfedgeBag h_bag(ctx.hex, selected_halfedges);
	auto selected_facets = h_bag.get_facets_indexes();	

	FacetBag f_bag(ctx.hex, selected_facets);
	auto selected_edges = f_bag.get_edges();

	const int n_layers = 2;

	int n_points = selected_edges.size();
	int n_cells = selected_facets.size();

	ctx.hex_bound->clear_surface();

	// Modify mesh
	int v_off = ctx.hex_bound->hex.points.create_points(n_points);
	int c_off = ctx.hex_bound->hex.create_cells(n_cells * n_layers);

	std::vector<std::vector<int>> polys(ctx.hex.nverts(), std::vector<int>(n_layers + 1, -1));

	for (int i = 0; i < static_cast<int>(selected_edges.size()); ++i) {
		auto e = Volume::Halfedge(ctx.hex, selected_edges[i]);
		ctx.hex_bound->hex.points[v_off + i] = (e.from().pos() + e.to().pos()) / 2.0;

		polys[e.from()][0] = e.from();
		polys[e.from()][1] = v_off + i;
		polys[e.from()][2] = e.to();
	}

	std::vector<int> count(ctx.hex.ncells(), false);

	auto selected_c_f = h_bag.get_cells_and_facets();
	auto &emb_attr = *ctx.emb_attr;
	for (int i = 0; i < n_cells; ++i) {

		Volume::Cell c(ctx.hex, selected_c_f[i][0]);
		Volume::Facet sf(ctx.hex, selected_c_f[i][1]);

		for (int lv = 0; lv < 4; ++lv) {
			int flc = 0;
			for (int lc = 0; lc < 8; ++lc) {
				if (c.vertex(lc) == polys[sf.vertex(lv)][0]) {
					flc = lc;
					break;
				}
			}
			int flc2 = 0;
			for (int lc = 0; lc < 8; ++lc) {
				if (c.vertex(lc) == polys[sf.vertex(lv)][2]) {
					flc2 = lc;
					break;
				}
			}
			// First hex
			ctx.hex_bound->hex.vert(c_off + i * n_layers, flc) = polys[sf.vertex(lv)][0];
			ctx.hex_bound->hex.vert(c_off + i * n_layers, flc2) = polys[sf.vertex(lv)][1];
			// Second hex
			ctx.hex_bound->hex.vert(c_off + i * n_layers + 1, flc) = polys[sf.vertex(lv)][1];
			ctx.hex_bound->hex.vert(c_off + i * n_layers + 1, flc2) = polys[sf.vertex(lv)][2];
		}

		// Update embedding
		Volume::Cell new_cell_0(ctx.hex, c_off + i * n_layers);
		Volume::Cell new_cell_1(ctx.hex, c_off + i * n_layers + 1);
		
		for (int lf = 0; lf < 6; ++lf) {
			emb_attr[new_cell_0.facet(lf)] = emb_attr[c.facet(lf)];
			emb_attr[new_cell_1.facet(lf)] = emb_attr[c.facet(lf)];
		}

		++count[c];
	}

	std::vector<bool> cells_to_kill(ctx.hex.ncells(), false);
	for (auto cf : selected_c_f) {
		cells_to_kill[cf[0]] = true;
	}

	ctx.hex.disconnect();
	ctx.hex.delete_cells(cells_to_kill);
	// Recompute mesh
	ctx.hex.connect();
	// Recompute layers
	layer_selector.init();

	// Update view
	ctx.recompute_hex();

	// Clear
	clear();
}

bool HexSplitTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void HexSplitTool::escape_callback() {

}