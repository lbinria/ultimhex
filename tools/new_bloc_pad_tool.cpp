#include "new_bloc_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"





bool LayerPad2::draw_object_properties() {
	if(ImGui::Button("Layer pad 2 !")) {
		ctx.gui_mode = NewBlocPadding;
		return true;
	}

	return false;
}

void LayerPad2::draw_viewer_properties() {}

void LayerPad2::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
	gl_draw::draw_path(selected_path, selected_color, false);

}

void LayerPad2::hover_callback(double x, double y, int source) {
	
	if (ctx.left_mouse_pressed && ctx.hex.connected() && ctx.is_cell_hovered() && ctx.is_cell_edge_hovered() && ctx.is_cell_lfacet_hovered()) {
		selected_path.clear();

		Volume::Cell um_c(ctx.hex, ctx.hovered_cell);
		Volume::Halfedge hovered_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.hovered_edge, ctx.hovered_cell_lfacet)));
		
		if (hovered_he.active()) {
			auto n_he = hovered_he.next();
			auto p_he = hovered_he.prev();

			auto n_dir = n_he.to().pos() - n_he.from().pos();
			auto p_dir = p_he.from().pos() - p_he.to().pos();

			auto a = hovered_he.from().pos() + p_dir * 0.2;
			auto b = hovered_he.to().pos() +  n_dir * 0.2;

			selected_path.push_back(a);
			selected_path.push_back(b);
		}

	}

}

void LayerPad2::mouse_button_callback(int button, int action, int mods, int source) {



}

void LayerPad2::scroll_callback(double xoffset, double yoffset) {

}



void LayerPad2::validate_callback() {

}

bool LayerPad2::is_compatible() { 
	return true;
}

void LayerPad2::escape_callback() {
	// Clear tool
	clear();
}

