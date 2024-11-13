#include "hover_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool HoverTool::draw_object_properties() {
	if(ImGui::Button("Hover")) {
		ctx.gui_mode = Hover;
		return true;
	}

	ImGui::Separator();

	return false;
}

void HoverTool::draw_viewer_properties() {
	ImGui::Checkbox("Show cell overlay", &ctx.show_hovered_cell_overlay_);
	ImGui::Checkbox("Show cell facet overlay", &ctx.show_hovered_cell_facet_overlay_);
	ImGui::SliderFloat("Thickness", &ctx.overlay_thickness, 1., 5.);
}

void HoverTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	if (ctx.show_hovered_cell_facet_overlay_) {
		if (ctx.is_cell_hovered() && ctx.is_cell_facet_hovered()) {
			gl_draw::draw_cell_facet_overlay(ctx.mesh_, ctx.hovered_cell, ctx.hovered_cell_lfacet, colorMapInfo, 0.0, ctx.overlay_thickness);	
		}
	}
	// Cell
	if (ctx.show_hovered_cell_overlay_) {
		if (ctx.is_cell_hovered()) {
			gl_draw::draw_cell_overlay(ctx.mesh_, ctx.hovered_cell, colorMapInfo, 0.0, ctx.overlay_thickness);
		}
	}
	if (ctx.is_cell_selected() && ctx.is_cell_lfacet_selected()) {
		gl_draw::draw_cell_facet_overlay(ctx.mesh_, ctx.selected_cell, ctx.selected_cell_lfacet, colorMapInfo, 0.5);
	}
	if (ctx.is_cell_selected()) {
		gl_draw::draw_cell_overlay(ctx.mesh_, ctx.selected_cell, colorMapInfo, 0.5);
	}
}

void HoverTool::hover_callback(double x, double y, int source) {

}

void HoverTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void HoverTool::scroll_callback(double xoffset, double yoffset) {}

void HoverTool::validate_callback() {

}

void HoverTool::escape_callback() {
	
}

bool HoverTool::is_compatible() { 
	return true;
}