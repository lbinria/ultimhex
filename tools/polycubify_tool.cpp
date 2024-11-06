#include "polycubify_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool PolycubifyTool::draw_object_properties() {
	int nhex_wanted = 3000;

	if (ImGui::Button("Polycubify !")) {

		// Get UM cell facet attribute tet_flag from GEO mesh
		UM::CellFacetAttribute<int> tet_flag(ctx.tet, -1);
		um_bindings::um_attr_from_geo_attr<GEO::MESH_CELL_FACETS>(ctx.mesh_, "tet_flag", ctx.tet, tet_flag.ptr);

		try {
			BenjaminAPI::polycubify(ctx.tet, tet_flag, ctx.hex, nhex_wanted);
		} catch (const std::runtime_error &e) {
			Logger::warn("An error occur when trying to polycubify. Detail: " + std::string(e.what()));
			std::cout << "polycubify fail" << std::endl;
			return false;
		}

		HexBoundary hex_bound(ctx.hex);
		// Replace current GEO mesh by UM Hex
		um_bindings::geo_mesh_from_hexboundary(hex_bound, ctx.mesh_);

		// TODO encapsulate in atomic unit ! + try catch to guarentee consistency

		// Write mesh

		// Save mesh metadata in json !!!!
		ctx.mesh_metadata = { 
			.filename = "polycubified.geogram", 
			.cell_type = GEO::MESH_HEX, 
			.attributes = {} 
		};
		write_by_extension(ctx.mesh_metadata.filename, hex_bound.hex, {{}, {}, {}, {}});
		ctx.mesh_metadata.save();

		// View
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		ctx.view.show_surface_ = false;
		ctx.view.show_volume_ = true;

		return true;
	}
	return false;
}

void PolycubifyTool::draw_viewer_properties() {

}

void PolycubifyTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void PolycubifyTool::hover_callback(double x, double y, int source) {

}

void PolycubifyTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void PolycubifyTool::validate_callback() {

}

void PolycubifyTool::escape_callback() {
	
}

bool PolycubifyTool::is_compatible() { 
	
	auto mesh_metadata_attr = ctx.mesh_metadata.get_attr("tet_flag");

	return ctx.mesh_metadata.cell_type == MESH_TET 
		&& mesh_metadata_attr.has_value() 
		&& mesh_metadata_attr.value().where == GEO::MESH_CELL_FACETS;
}