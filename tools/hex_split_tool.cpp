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
		
		ctx.view.attribute_subelements_ = GEO::MeshElementsFlags::MESH_NONE;
		ctx.view.cells_shrink_ = 0.0f;
		ctx.view.switch_to_volume_select_mode();
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
	
	const auto &s_halfedges = layer_selector.get_selected_halfedges();

	std::vector<Volume::Halfedge> selected_halfedges;
	for (auto hi : s_halfedges) {
		selected_halfedges.push_back(Volume::Halfedge(ctx.hex_bound->hex, hi));
	}

	HalfedgeBag h_bag(selected_halfedges);
	auto selected_facets = h_bag.get_facets_indexes();
	auto selected_facets2 = h_bag.get_facets_indexes2();

	FacetBag f_bag(ctx.hex, selected_facets);
	auto selected_edges = f_bag.get_edges();

	int n_layers = 1;

	int n_points = selected_edges.size();
	int n_cells = selected_facets.size();

	ctx.hex_bound->clear_surface();

	// Modify mesh
	int v_off = ctx.hex_bound->hex.points.create_points(n_points);
	int c_off = ctx.hex_bound->hex.create_cells(n_cells * (n_layers + 1));

	std::vector<int> edge_2_pts_indexes(ctx.hex.nverts(), -1);

	for (int i = 0; i < static_cast<int>(selected_edges.size()); ++i) {
		auto e = Volume::Halfedge(ctx.hex, selected_edges[i]);
		ctx.hex_bound->hex.points[v_off + i] = (e.from().pos() + e.to().pos()) / 2.0;
		edge_2_pts_indexes[e.from()] = v_off + i;
	}


	int llv[4] = {0,1,3,2};
	int llv2[4] = {2,3,1,0};

	for (int ci = 0; ci < static_cast<int>(selected_facets2.size()); ++ci) {
		auto start_f = Volume::Facet(ctx.hex, selected_facets2[ci].first);
		auto end_f = Volume::Facet(ctx.hex, selected_facets2[ci].second);
		
		for (int lv = 0; lv < 4; ++lv) {
			ctx.hex_bound->hex.vert(c_off + ci * 2, llv[lv]) = start_f.halfedge(lv).from();
		}

		for (int lv = 0; lv < 4; ++lv) {
			auto h = start_f.halfedge(lv);
			auto oh = h.opposite_f().next();
			ctx.hex_bound->hex.vert(c_off + ci * 2, llv[lv] + 4) = edge_2_pts_indexes[oh.from()];
			ctx.hex_bound->hex.vert(c_off + ci * 2 + 1, llv[lv]) = edge_2_pts_indexes[oh.from()];

		}

		for (int lv = 0; lv < 4; ++lv) {
			ctx.hex_bound->hex.vert(c_off + ci * 2 + 1, llv2[lv] + 4) = end_f.halfedge(lv).to();
		}
	}


	std::vector<bool> cells_to_kill(ctx.hex.ncells(), false);
	auto selected_cells = h_bag.get_cells_indexes();
	for (auto c : selected_cells) {
		cells_to_kill[c] = true;
	}

	// Update embedding
	auto &emb_attr = *ctx.emb_attr;
	{
		bool done = false;
		while (!done) {
			done = true;
			for (auto f : ctx.hex.iter_facets()) {
				// If cell is to kill and facet has an embedding (>= 0)
				if (emb_attr[f] < 0 || !cells_to_kill[f.cell()]) 
					continue;
				
				// We move its embedding to the opposite interior facet
				auto f_next = f.halfedge(0).opposite_f().next().next().opposite_f().facet().opposite();
				
				if (f_next.active()) {
					std::swap(emb_attr[f], emb_attr[f_next]);
					done = false;
				}
			}
		}
	}


	ctx.hex.disconnect();
	ctx.hex.delete_cells(cells_to_kill);
	ctx.hex.delete_isolated_vertices();
	// Recompute mesh
	ctx.hex.connect();
	// Recompute layers
	layer_selector.init();

	// Update view
	ctx.recompute_hex();

	// Clear
	clear();


	// write_by_extension("hex_split.geogram", ctx.hex_bound->hex, {});

}

bool HexSplitTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void HexSplitTool::escape_callback() {

}