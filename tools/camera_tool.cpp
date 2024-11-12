#include "camera_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool CameraTool::draw_object_properties() {
	if(ImGui::Button("Camera")) {
		ctx.gui_mode = Camera;
		return true;
	}

	ImGui::Separator();

	return false;
}

void CameraTool::draw_viewer_properties() {

}

void CameraTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void CameraTool::hover_callback(double x, double y, int source) {

}

void CameraTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void CameraTool::scroll_callback(double xoffset, double yoffset) {}

void CameraTool::validate_callback() {

}

void CameraTool::escape_callback() {
	
}

bool CameraTool::is_compatible() { 
	return true;
}