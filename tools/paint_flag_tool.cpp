#include "paint_flag_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../tag_face.h"


bool PaintFlagTool::draw_object_properties() {

	// auto t = convert_to_ImTextureID(colormaps_[current_colormap_index_].texture);
	// ImGui::ColorButton("+X", ImVec4(1,0,0,1));

	ImGui::TextUnformatted("Paint flags");

	if(ImGui::Button("No Paint")) {
		paint_value = -1;
		ctx.gui_mode = Painting;
	}

	if(ImGui::Button("-X")) {
		paint_value = 0;
		ctx.gui_mode = Painting;
	}
	ImGui::SameLine();
	if(ImGui::Button("-Y")) {
		paint_value = 1;
		ctx.gui_mode = Painting;
	}
	ImGui::SameLine();
	if(ImGui::Button("-Z")) {
		paint_value = 2;
		ctx.gui_mode = Painting;
	}

	if(ImGui::Button("+X")) {
		paint_value = 3;
		ctx.gui_mode = Painting;
	}
	ImGui::SameLine();
	if(ImGui::Button("+Y")) {
		paint_value = 4;
		ctx.gui_mode = Painting;
	}
	ImGui::SameLine();
	if(ImGui::Button("+Z")) {
		paint_value = 5;
		ctx.gui_mode = Painting;
	}	

	ImGui::Separator();

	if (ImGui::Button("Compute flags !")) {
		
		// Compute flag on tet and tri
		// TetBoundary tet_bound(tet);
		ctx.tet_bound.update();
		
		// To GEO mesh
		// TODO maybe move under tet_bound.set_attribute_to_surface(tet_flag, tri_flag);
		um_bindings::geo_mesh_from_tetboundary(ctx.tet_bound, ctx.mesh_);

		UM::CellFacetAttribute<int> tet_flag(ctx.tet, -1);
		UM::FacetAttribute<int> tri_flag(ctx.tet_bound.tri, -1);

		// Compute flag
		algo::naive_tag(ctx.tet, tet_flag);
		// Transfert flag from tet to tri for display
		ctx.tet_bound.set_attribute_to_surface(tet_flag, tri_flag);
		// Update GEO mesh attribute "flag"
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(ctx.tet, tet_flag.ptr, "tet_flag", ctx.mesh_);
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.tet_bound.tri, tri_flag.ptr, "flag", ctx.mesh_);


		// TODO encapsulate in atomic unit ! + try catch to guarentee consistency
		// Save mesh metadata
		ctx.mesh_metadata = { 
			.filename = "flagged.geogram", 
			.cell_type = GEO::MESH_TET, 
			.attributes = {
				{
					.name = "tet_flag", 
					.type = "int", 
					.where = MESH_CELL_FACETS
				}
			} 
		};
		// Write mesh
		write_by_extension(ctx.mesh_metadata.filename, ctx.tet_bound.tet, {{}, {}, {{"tet_flag", tet_flag.ptr}}, {}});
		ctx.mesh_metadata.save();

		ctx.view.show_surface_ = true;
		ctx.view.show_volume_ = false;
		
	}

	// if (ImGui::Button("Compute patches !")) {
	// 	// Compute flag on tet and tri
	// 	TetBoundary tet_bound(ctx.tet);
	// 	UM::FacetAttribute<int> tri_flag(tet_bound.tri, -1);
	// 	um_bindings::um_attr_from_geo_attr<GEO::MESH_FACETS>(mesh_, "flag", tet_bound.tri, tri_flag.ptr);
	// 	flag_dirs = compute_patches(tet_bound.tri, tri_flag);
	// }

	ImGui::Separator();
}

void PaintFlagTool::draw_viewer_properties() {}

void PaintFlagTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	// Nothing to draw
}

void PaintFlagTool::hover_callback(double x, double y, int source) {

	std::cout << "left mouse : " << ctx.left_mouse_pressed << std::endl;
	// Just when left mouse is currently pressed (drag)
	if (!(ctx.left_mouse_pressed && ctx.is_cell_facet_hovered()))
		return;

	GEO::Attribute<GEO::signed_index_t> flag(
		ctx.mesh_.facets.attributes(), "flag"
	);

	flag[ctx.hovered_facet] = paint_value;
}

void PaintFlagTool::mouse_button_callback(int button, int action, int mods, int source) {
	// Transfert attribute from surface tri to volume tet
	FacetAttribute<int> tri_flag(ctx.tet_bound.tri, -1);
	CellFacetAttribute<int> tet_flag(ctx.tet_bound.tet, -1);
	um_bindings::um_attr_from_geo_attr<GEO::MESH_FACETS>(ctx.mesh_, "flag", ctx.tet_bound.tri, tri_flag.ptr);

	ctx.tet_bound.set_attribute_to_volume(tri_flag, tet_flag);
	um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(ctx.tet_bound.tet, tet_flag.ptr, "tet_flag", ctx.mesh_);
}

void PaintFlagTool::scroll_callback(double xoffset, double yoffset) {}

void PaintFlagTool::validate_callback() {

}

bool PaintFlagTool::is_compatible() { 
	return !ctx.mesh_metadata.filename.empty() && ctx.mesh_metadata.cell_type == MESH_TET;
}

void PaintFlagTool::escape_callback() {

}

