#include "filter_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

std::vector<int> myloop(UM::Hexahedra &hex, UM::Volume::Facet &f) {
	assert(hex.connected());

	std::vector<int> facets;

	std::vector<bool> visited(hex.nfacets(), false);
	std::queue<int> q;
	
	q.push(f);
	visited[f] = true;

	while (!q.empty()) {

		auto f_idx = q.front();
		auto f = Volume::Facet(hex, f_idx);
		q.pop();

		for (int i = 0; i < f.nhalfedges(); i++) {
			auto h = f.halfedge(i);

			auto opp_c = h.next().opposite_f().opposite_c();
			if (!opp_c.active())
				continue;

			auto nxt_f = opp_c.opposite_f().facet();
			if (!visited[nxt_f]) {
				q.push(nxt_f);
				visited[nxt_f] = true;
			}
		}

		facets.push_back(f);
		
	}

	return facets;
}

void FilterTool::reset_filters() {
	GEO::Attribute<bool> facet_filter(
		ctx.mesh_.facets.attributes(), "filter"
	);

	for (auto f : ctx.mesh_.facets)
		facet_filter[f] = true;

	GEO::Attribute<bool> cell_filter(
		ctx.mesh_.cells.attributes(), "filter"
	);

	for (auto c : ctx.mesh_.cells) 
		cell_filter[c] = true;

	// ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_FACETS);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
}

bool FilterTool::draw_object_properties() {
		


		if (ImGui::Button("Reset filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
			reset_filters();
		}

		if (ImGui::Button("Switch filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {

			GEO::Attribute<bool> filter(
				ctx.mesh_.facets.attributes(), "filter"
			);

			for (auto f : ctx.mesh_.facets)
				filter[f] = !filter[f];

			ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_FACETS);
		}

		if (ImGui::Button("Chart filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
			reset_filters();
			ctx.gui_mode = Filter;
			mode = Chart;
		}

		if (ImGui::Button("Hex minecraft filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
			reset_filters();
			ctx.gui_mode = Filter;
			mode = Minecraft;
		}

		if (ImGui::Button("Hex layer filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
			reset_filters();
			ctx.gui_mode = Filter;
			mode = Layer;
		}

		if (mode == Minecraft) {

			// Paint mode combo box selection
			ImGui::Text("Mode");
			if (ImGui::BeginCombo("a", minecraft_modes[current_minecraft_mode])) {

				for (int i = 0; i < IM_ARRAYSIZE(minecraft_modes); i++)
				{
					bool isSelected = (current_minecraft_mode == i);
					if (ImGui::Selectable(minecraft_modes[i], isSelected))
					{
						current_minecraft_mode = i;
					}

					if (isSelected)
						ImGui::SetItemDefaultFocus();
				}

				ImGui::EndCombo();
			}
		}

	// ImVec2 window_pos = ImGui::GetWindowPos(); 
	// ImVec2 window_size = ImGui::GetWindowSize(); 
	ImVec2 mouse_pos = ImGui::GetMousePos();
	// ImVec2 cursor_pos = ImVec2(window_pos.x + mouse_pos.x, window_pos.y + mouse_pos.y); 
	// ImGui::GetBackgroundDrawList()->AddCircle(cursor_pos, window_size.x * 0.6f, IM_COL32(255, 0, 0, 200), 0, 5); 
	// ImGui::GetBackgroundDrawList()->AddCircle(mouse_pos, ctx.brush_size, IM_COL32(230, 180, 255, 250), 0, 3); 
	ImGui::GetBackgroundDrawList()->AddRect(ImVec2(mouse_pos.x - ctx.brush_size / 2, mouse_pos.y - ctx.brush_size / 2), ImVec2(mouse_pos.x + ctx.brush_size / 2, mouse_pos.y + ctx.brush_size / 2), IM_COL32(230, 180, 255, 250), 0.2, 0, 3); 
	// ImGui::GetForegroundDrawList()->AddCircle(window_center, window_size.y * 0.6f, IM_COL32(0, 255, 0, 200), 0, 10);

	return false;
}

void FilterTool::draw_viewer_properties() {



}

void FilterTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
	for (auto c : selected_cells)
		gl_draw::draw_cell_overlay(ctx.mesh_, c, colorMapInfo, 0.5);

	for (auto c : test)
		gl_draw::draw_cell_overlay(ctx.mesh_, c, colorMapInfo, 0.5);
}

void FilterTool::hover_callback(double x, double y, int source) {

	ctx.region_selection_activated = ctx.left_mouse_pressed;

	if (ctx.left_mouse_pressed && ctx.is_cell_hovered() && mode == Minecraft) {
		minecraft();
	}

}

void FilterTool::mouse_button_callback(int button, int action, int mods, int source) {
	
	if (ctx.is_facet_hovered() && mode == Chart) {
		filter_chart();

	} else if (mode == Minecraft) {

		// Filter all selected cells
		GEO::Attribute<bool> cell_filter(
			ctx.mesh_.cells.attributes(), "filter"
		);

		if (current_minecraft_mode == Dig) {
			// Just filter all selected cells
			for (auto c : selected_cells)
				cell_filter[c] = false;

			selected_cells.clear();

		} else if (current_minecraft_mode == Levelling) {
			// We take first cell as levelling reference
			// Then we search for 

		}


		ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_CELLS);

	} else if (mode == Layer && ctx.is_cell_facet_hovered()) {

		// Filter all selected cells
		GEO::Attribute<bool> cell_filter(
			ctx.mesh_.cells.attributes(), "filter"
		);

		test.clear();
		auto f = Volume::Facet(ctx.hex, ctx.um_hovered_cell_facet());
		auto facets = myloop(ctx.hex, f);
		for (auto f : facets) {
			Volume::Facet cur_f(ctx.hex, f);

			cell_filter[cur_f.cell()] = false;
			// test.push_back(cur_f.cell());
		}


		ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_CELLS);

		// // Filter all selected cells
		// GEO::Attribute<bool> cell_filter(
		// 	ctx.mesh_.cells.attributes(), "filter"
		// );

		// Volume::Facet um_hovered_f(ctx.hex, ctx.um_hovered_cell_facet());
		// auto start_h = um_hovered_f.halfedge(0).opposite_f().next();
		// loop_cutyo(ctx.hex, start_h, [&](Volume::Facet &f) {
		// 	// um_bindings::geo_cell_index_from_facet_index(f);
		// 	cell_filter[f.cell()] = false;
		// });

		// ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_CELLS);

	}

}

void FilterTool::minecraft() {

	selected_cells.insert(ctx.hovered_cells.begin(), ctx.hovered_cells.end());

}

void FilterTool::filter_chart() {


	std::vector<int> facet_by_features;

	// TODO duplicate code in paint flag tool
	std::vector<bool> is_feature(ctx.tet_bound->tri.ncorners(), false);

	// Group halfedge by vertices
	std::map<std::pair<int, int>, int> halfedge_by_vertices;
	for (auto h : ctx.tet_bound->tri.iter_halfedges()) {
		halfedge_by_vertices[{h.from(), h.to()}] = h;
	}

	for (auto e : ctx.mesh_.edges) {
		auto v0 = ctx.mesh_.edges.vertex(e, 0);
		auto v1 = ctx.mesh_.edges.vertex(e, 1);

		int h1 = halfedge_by_vertices[{v0, v1}];
		int h2 = halfedge_by_vertices[{v1, v0}];
		is_feature[h1] = true;
		is_feature[h2] = true;
	}

	DisjointSet ds(ctx.tet_bound->tri.nfacets());
	for (auto h : ctx.tet_bound->tri.iter_halfedges()) {
		auto opp_h = h.opposite();

		if (!opp_h.active() || is_feature[opp_h])
			continue;

		ds.merge(h.facet(), opp_h.facet());
	}
	facet_by_features.clear();
	facet_by_features.resize(ctx.tet_bound->tri.nfacets());
	ds.get_sets_id(facet_by_features);

	int feature = facet_by_features[ctx.selected_facet];

	GEO::Attribute<bool> filter(
		ctx.mesh_.facets.attributes(), "filter"
	);

	for (int i = 0; i < facet_by_features.size(); i++) {
		filter[i] = facet_by_features[i] != feature;
	}

	ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_FACETS);
}

void FilterTool::scroll_callback(double xoffset, double yoffset) {
	if (yoffset != 0)
		ctx.brush_size += (int)yoffset;
}

void FilterTool::validate_callback() {

}

void FilterTool::escape_callback() {
	
}

bool FilterTool::is_compatible() { 
	return true;
}