#include "embedit_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool EmbeditTool::draw_object_properties() {


	// if(ImGui::Button("Embedit##btn_embedit_embedit_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {


	// 	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);

	// 	FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
	// 	ctx.hex_bound->set_attribute_to_surface(*ctx.emb_attr, surf_emb_attr);
	// 	um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.hex_bound->quad, surf_emb_attr.ptr, "emb", ctx.mesh_);

	// 	// Search for min / max value of embedding attribute
	// 	int min = std::numeric_limits<int>::max(), max = 0;
	// 	for (auto f : ctx.hex_bound->quad.iter_facets()) {
	// 		if (surf_emb_attr[f] < min)
	// 			min = surf_emb_attr[f];
	// 		if (surf_emb_attr[f] > max)
	// 			max = surf_emb_attr[f];
	// 	}

	// 	// Display surface with embedding attr
	// 	ctx.view.change_mode(ViewBinding::Mode::Surface);
	// 	ctx.view.attribute_subelements_ = GEO::MeshElementsFlags::MESH_FACETS;
	// 	ctx.view.attribute_ = "facets.emb";
	// 	ctx.view.attribute_name_ = "emb";
	// 	ctx.view.attribute_min_ = min;
	// 	ctx.view.attribute_max_ = max;

	// 	ctx.gui_mode = GUIMode::Embedit;

	// 	return true;
	// }

	if (!is_init && ImGui::Button("Init embedit##btn_charts_embedit_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {

		
		// auto &tri_chart = *ctx.tri_chart;
		auto &quad_chart = *ctx.quad_chart;

		// um_bindings::geo_mesh_from_tetboundary(*ctx.tet_bound, ctx.mesh_);
		// um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.tet_bound->tri, tri_chart.ptr, "charts", ctx.mesh_);

		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.hex_bound->quad, quad_chart.ptr, "charts", ctx.mesh_);

		// // Search for min / max value of charts attribute
		// int min = std::numeric_limits<int>::max(), max = 0;
		// for (auto f : ctx.tet_bound->tri.iter_facets()) {
		// 	if (tri_chart[f] < min)
		// 		min = tri_chart[f];
		// 	if (tri_chart[f] > max)
		// 		max = tri_chart[f];
		// }

		// Search for min / max value of charts attribute
		int min = std::numeric_limits<int>::max(), max = 0;
		for (auto f : ctx.hex_bound->quad.iter_facets()) {
			if (quad_chart[f] < min)
				min = quad_chart[f];
			if (quad_chart[f] > max)
				max = quad_chart[f];
		}

		// Number of charts is max - min, with min included => + 1 
		n_charts = max - min + 1;
		std::cout << "n charts: " << n_charts << std::endl;

		// Display surface with embedding attr
		ctx.view.change_mode(ViewBinding::Mode::Surface);
		ctx.view.attribute_subelements_ = GEO::MeshElementsFlags::MESH_FACETS;
		ctx.view.attribute_ = "facets.charts";
		ctx.view.attribute_name_ = "charts";
		ctx.view.attribute_min_ = min;
		ctx.view.attribute_max_ = max;

		ctx.gui_mode = GUIMode::Embedit;
		mode = Paint;

		is_init = true;
		return true;
	}

	if (!is_init)
		return false;

	if (ImGui::Button("Pick##btn_charts_pick_embedit_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		mode = Pick;
	}
	if (ImGui::Button("Paint##btn_charts_paint_embedit_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		mode = Paint;
	}


	if (ImGui::Button("Apply##btn_charts_apply_embedit_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		BenjaminAPI::embeditapply(ctx.hex_bound->hex, *ctx.emb_attr, ctx.hex_bound->quad, *ctx.quad_chart, ctx.tet_bound->tri, *ctx.tri_chart);
		// BenjaminAPI::smooth(ctx.hex_bound->hex, *ctx.emb_attr, ctx.tet_bound->tri, *ctx.tri_chart);

		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.view.change_mode(ViewBinding::Mode::Volume);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);


	}

	// for (int i = 0; i < n_charts; i++) {
	// 	ImGui::PushID(i);
	// 	auto btn_color = ctx.view.getColorFromTexture1D(ctx.view.colormaps[ctx.view.current_colormap_index_], i);
	// 	ImGui::PushStyleColor(ImGuiCol_Button, btn_color);
	// 	if (ImGui::Button(" ", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
	// 		selected_chart_color = i;
	// 	}

	// 	ImGui::PopStyleColor();
	// 	ImGui::PopID();
	// }

	return false;
}

void EmbeditTool::draw_viewer_properties() {

}

void EmbeditTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void EmbeditTool::hover_callback(double x, double y, int source) {

	// Paint chart color
	if (mode == Paint && ctx.left_mouse_pressed && ctx.is_facet_hovered()) {

		GEO::Attribute<GEO::signed_index_t> charts(
			ctx.mesh_.facets.attributes(), "charts"
		);
		charts[ctx.hovered_facet] = selected_chart_color;
		// (*ctx.tri_chart)[ctx.hovered_facet] = selected_chart_color;
		(*ctx.quad_chart)[ctx.hovered_facet] = selected_chart_color;
	}

}

void EmbeditTool::mouse_button_callback(int button, int action, int mods, int source) {
	if (ctx.is_facet_hovered() && mode == Pick) {
		selected_chart_color =  (*ctx.quad_chart)[ctx.hovered_facet];
	}
}

void EmbeditTool::scroll_callback(double xoffset, double yoffset) {}

void EmbeditTool::validate_callback() {

}

void EmbeditTool::escape_callback() {
	
}

bool EmbeditTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}