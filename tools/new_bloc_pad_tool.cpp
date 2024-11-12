#include "new_bloc_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"





bool NewBlocPadTool::draw_object_properties() {
	if(ImGui::Button("New Bloc padding")) {
		ctx.gui_mode = NewBlocPadding;
		return true;
	}

	return false;
}

void NewBlocPadTool::draw_viewer_properties() {}

void NewBlocPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
	// glupDisable(GLUP_LIGHTING);
	// glupMatrixMode(GLUP_PROJECTION_MATRIX);
	// glupLoadIdentity();
	// glupOrtho2D(-100,100,-100,100);

	// glupMatrixMode(GLUP_MODELVIEW_MATRIX);
	// glupLoadIdentity();
	
	
	// // glupPushMatrix();

	// glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(1.f, 0.f,0.f, 0.5f).data());
	// glupBegin(GLUP_QUADS);
	// glupPrivateVertex2f(0,0);
	// glupPrivateVertex2f(20,0);
	// glupPrivateVertex2f(20,20);
	// glupPrivateVertex2f(0,20);

	// glupEnd();

	// // glupPopMatrix();
	// glupEnable(GLUP_LIGHTING);

}

void NewBlocPadTool::hover_callback(double x, double y, int source) {

}

void NewBlocPadTool::mouse_button_callback(int button, int action, int mods, int source) {
	
}

void NewBlocPadTool::scroll_callback(double xoffset, double yoffset) {

}



void NewBlocPadTool::validate_callback() {

}

bool NewBlocPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void NewBlocPadTool::escape_callback() {
	// Clear tool
	clear();
}

