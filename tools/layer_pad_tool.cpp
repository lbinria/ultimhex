#include "layer_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

#include "../helpers.h"

bool LayerPadTool::draw_object_properties() {
	if(ImGui::Button("Init##layer_pad_tool_init")) {

		ctx.gui_mode = LayerPadding;		
		return true;
	}

	return false;
}

void LayerPadTool::draw_viewer_properties() {}

void LayerPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void LayerPadTool::hover_callback(double x, double y, int source) {

}

void LayerPadTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void LayerPadTool::scroll_callback(double xoffset, double yoffset) {}

void LayerPadTool::validate_callback() {

}

bool LayerPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}