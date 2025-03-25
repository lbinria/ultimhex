#include "smooth_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool SmoothTool::draw_object_properties() {

	if(ImGui::Button("Smooth##btn_smooth_smooth_tool", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		
		// Get alias
		CellFacetAttribute<int> &emb_attr = *ctx.emb_attr;

		// Check validity of emb
		for (auto f : ctx.hex_bound->hex.iter_facets()) if (f.on_boundary() && emb_attr[f] < 0) {
			for (auto in_f : ctx.hex_bound->hex.iter_facets()) if (emb_attr[in_f] >= 0 || !in_f.on_boundary()) emb_attr[in_f] = 0;
			std::cout << "bad embedding !" << std::endl;
			// TODO should visualize feedback !!!
			FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
			ctx.hex_bound->set_attribute_to_surface(emb_attr, surf_emb_attr);
			write_by_extension("bad_emb.geogram", ctx.hex_bound->quad, {{}, {{"emb", surf_emb_attr.ptr}}, {}});
			
			return false;
		}


		BenjaminAPI::smooth(ctx.hex_bound->hex, emb_attr, ctx.tet_bound->tri, *ctx.tri_chart);

		ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

		ctx.gui_mode = GUIMode::Smooth;

		return true;
	}

	return false;
}

void SmoothTool::draw_viewer_properties() {

}

void SmoothTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void SmoothTool::hover_callback(double x, double y, int source) {

}

void SmoothTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void SmoothTool::scroll_callback(double xoffset, double yoffset) {}

void SmoothTool::validate_callback() {

}

void SmoothTool::escape_callback() {
	
}

bool SmoothTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}