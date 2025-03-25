#include "embedit_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool EmbeditTool::draw_object_properties() {


	if(ImGui::Button("Embedit##btn_embedit_embedit_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {

		// Transfert attribute from UM to GEO
		// um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.hex_bound->quad, ctx.emb_attr.get()->ptr, "emb", ctx.mesh_);

		// TODO here should transfert emb from hex cellfacet to quad facet...
		FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
		ctx.hex_bound->set_attribute_to_surface(*ctx.emb_attr, surf_emb_attr);
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.hex_bound->quad, surf_emb_attr.ptr, "emb", ctx.mesh_);

		// Display surface with embedding attr
		ctx.view.change_mode(ViewBinding::Mode::Surface);
		ctx.view.attribute_ = "facets.emb";
		ctx.view.attribute_name_ = "emb";
		// ctx.view.attribute_min_ = 0;
		// ctx.view.attribute_max_ = 2;

		ctx.gui_mode = GUIMode::Embedit;

		return true;
	}

	return false;
}

void EmbeditTool::draw_viewer_properties() {

}

void EmbeditTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void EmbeditTool::hover_callback(double x, double y, int source) {

}

void EmbeditTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void EmbeditTool::scroll_callback(double xoffset, double yoffset) {}

void EmbeditTool::validate_callback() {

}

void EmbeditTool::escape_callback() {
	
}

bool EmbeditTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}