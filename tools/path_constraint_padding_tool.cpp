#include "path_constraint_padding_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool PathConstraintPaddingTool::draw_object_properties() {
	if(ImGui::Button("Init", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		ctx.gui_mode = Hover;
		return true;
	}

	ImGui::Separator();

	return false;
}

void PathConstraintPaddingTool::draw_viewer_properties() {

}

void PathConstraintPaddingTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void PathConstraintPaddingTool::hover_callback(double x, double y, int source) {

}

void PathConstraintPaddingTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void PathConstraintPaddingTool::scroll_callback(double xoffset, double yoffset) {}

void PathConstraintPaddingTool::validate_callback() {

}

void PathConstraintPaddingTool::escape_callback() {
	
}

bool PathConstraintPaddingTool::is_compatible() { 
	return true;
}