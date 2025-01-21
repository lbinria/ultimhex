#include "paint_flag_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../tag_face.h"


void PaintFlagTool::init() {
	// Compute flag on tet and tri
	ctx.tet_bound.update();

	UM::CellFacetAttribute<int> tet_flag(ctx.tet, -1);
	UM::FacetAttribute<int> tri_flag(ctx.tet_bound.tri, -1);

	// Transfert flag from tet to tri for display
	ctx.tet_bound.set_attribute_to_surface(tet_flag, tri_flag);

	// Keep features
	std::vector<std::pair<int, int>> edges(ctx.mesh_.edges.nb());
	for (int e = 0; e < ctx.mesh_.edges.nb(); e++) {
		auto v0 = ctx.mesh_.edges.vertex(e, 0);
		auto v1 = ctx.mesh_.edges.vertex(e, 1);
		auto p0 = ctx.mesh_.vertices.point(v0);
		auto p1 = ctx.mesh_.vertices.point(v1);
		edges.push_back({v0, v1});
	}

	// To GEO mesh (we have to do that to have surface that match to volume, thanks to TetBound)
	um_bindings::geo_mesh_from_tetboundary(ctx.tet_bound, ctx.mesh_);

	// Restore features
	ctx.mesh_.edges.create_edges(edges.size());
	for (int e = 0; e < edges.size(); e++) {
		ctx.mesh_.edges.set_vertex(e, 0, edges[e].first);
		ctx.mesh_.edges.set_vertex(e, 1, edges[e].second);
	}

	// Update GEO mesh attribute "flag"
	um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(ctx.tet, tet_flag.ptr, "tet_flag", ctx.mesh_);
	um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.tet_bound.tri, tri_flag.ptr, "flag", ctx.mesh_);

	ctx.view.change_mode(ViewBinding::Mode::Surface);
}

bool PaintFlagTool::draw_object_properties() {

	// auto t = convert_to_ImTextureID(colormaps_[current_colormap_index_].texture);
	// ImGui::ColorButton("+X", ImVec4(1,0,0,1));

	// if (ImGui::CollapsingHeader("Paint flag")) {

		if (is_init && ctx.gui_mode != Painting && ctx.switch_mode != Painting && ImGui::Button("Paint !")) {
			ctx.gui_mode = Painting;
		}

		if(!is_init && ImGui::Button("Init flag painting")) {

			init();

			// Init feature bucket

			std::vector<bool> is_feature(ctx.tet_bound.tri.ncorners(), false);

			// Group halfedge by vertices
			std::map<std::pair<int, int>, int> halfedge_by_vertices;
			for (auto h : ctx.tet_bound.tri.iter_halfedges()) {
				halfedge_by_vertices[{h.from(), h.to()}] = h;
			}



			for (auto e : ctx.mesh_.edges) {
				auto v0 = ctx.mesh_.edges.vertex(e, 0);
				auto v1 = ctx.mesh_.edges.vertex(e, 1);

				int h1 = halfedge_by_vertices[{v0, v1}];
				int h2 = halfedge_by_vertices[{v1, v0}];
				is_feature[h1] = true;
				is_feature[h2] = true;
			}
			
			DisjointSet ds(ctx.tet_bound.tri.nfacets());
			for (auto h : ctx.tet_bound.tri.iter_halfedges()) {
				auto opp_h = h.opposite();

				if (!opp_h.active() || is_feature[opp_h])
					continue;

				ds.merge(h.facet(), opp_h.facet());
			}
			facet_by_features.clear();
			facet_by_features.resize(ctx.tet_bound.tri.nfacets());
			ds.get_sets_id(facet_by_features);


			// TODO maybe move that : Push attribute to metadata
			ctx.mesh_metadata.attributes.push_back(MeshMetadata::MetadataAttribute{ .name = "tet_flag", .type = "int", .where = MESH_CELL_FACETS });

			ctx.gui_mode = Painting;
			is_init = true;
		}

		// Only display paint menu when painting was initialized
		if (!is_init || ctx.gui_mode != Painting && ctx.switch_mode != Painting)
			return false;

		// Remove all flags of current mesh
		if (ImGui::Button("Remove all flags")) {
			UM::CellFacetAttribute<int> tet_flag(ctx.tet, -1);
			UM::FacetAttribute<int> tri_flag(ctx.tet_bound.tri, -1);
			// Update GEO mesh attribute "flag"
			um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(ctx.tet, tet_flag.ptr, "tet_flag", ctx.mesh_);
			um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.tet_bound.tri, tri_flag.ptr, "flag", ctx.mesh_);

			ctx.view.change_mode(ViewBinding::Mode::Surface);
		}

		// Compute all flags of current mesh in naive mode
		if (ImGui::Button("Compute all flags")) {

			UM::CellFacetAttribute<int> tet_flag(ctx.tet, -1);
			UM::FacetAttribute<int> tri_flag(ctx.tet_bound.tri, -1);

			// Compute flag
			algo::naive_tag(ctx.tet, tet_flag);
			// algo::generate_naive_flagging(ctx.tet, tet_flag);
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

			ctx.view.change_mode(ViewBinding::Mode::Surface);
		}

		if (ImGui::Button("List not flagged")) {
			GEO::Attribute<GEO::signed_index_t> flag(
				ctx.mesh_.facets.attributes(), "flag"
			);

			for (auto f : ctx.mesh_.facets) {
				if (flag[f] < 0)
					std::cout << "facette " << f << "isn't flagged" << std::endl;
			}
		}

		// Paint mode combo box selection
		ImGui::Text("Mode");
		if (ImGui::BeginCombo("a", modes[current_mode])) {

			for (int i = 0; i < IM_ARRAYSIZE(modes); i++)
			{
				bool isSelected = (current_mode == i);
				if (ImGui::Selectable(modes[i], isSelected))
				{
					current_mode = i;
				}

				if (isSelected)
					ImGui::SetItemDefaultFocus();
			}

			ImGui::EndCombo();
		}

		// Paint algo combo box selection
		ImGui::Text("Algorithm");
		if (ImGui::BeginCombo("b", algos[current_algo])) {

			for (int i = 0; i < IM_ARRAYSIZE(algos); i++)
			{
				bool isSelected = (current_algo == i);
				if (ImGui::Selectable(algos[i], isSelected))
				{
					current_algo = i;
				}

				if (isSelected)
					ImGui::SetItemDefaultFocus();
			}

			ImGui::EndCombo();
		}

		ImGui::Separator();
		ImGui::Text("Algorithm options");

		if (current_algo == None) {
			if(ImGui::Button("No Paint")) {
				value = -1;
			}

			if(ImGui::Button("-X")) {
				value = 0;
			}
			ImGui::SameLine();
			if(ImGui::Button("-Y")) {
				value = 1;
			}
			ImGui::SameLine();
			if(ImGui::Button("-Z")) {
				value = 2;
			}

			if(ImGui::Button("+X")) {
				value = 3;
			}
			ImGui::SameLine();
			if(ImGui::Button("+Y")) {
				value = 4;
			}
			ImGui::SameLine();
			if(ImGui::Button("+Z")) {
				value = 5;
			}
		}

		if (current_algo == Naive) {
			ImGui::TextUnformatted("Components used for orientation detection:");
			ImGui::Checkbox("X", &naive_constraints[0]);
			ImGui::Checkbox("Y", &naive_constraints[1]);
			ImGui::Checkbox("Z", &naive_constraints[2]);
		}

		ImGui::Separator();

	ImGui::Separator();

	return true;
}

void PaintFlagTool::draw_viewer_properties() {}

void PaintFlagTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	// Nothing to draw
}

void PaintFlagTool::hover_callback(double x, double y, int source) {

	// Just when left mouse is currently pressed (drag)
	if (!ctx.is_facet_hovered())
		return;

	if (ctx.left_mouse_pressed && current_mode == Facet) {

		if (current_algo == None) {
			
			GEO::Attribute<GEO::signed_index_t> flag(
				ctx.mesh_.facets.attributes(), "flag"
			);

			flag[ctx.hovered_facet] = value;

		} else if (current_algo == Naive) {

			UM::Surface::Facet f(ctx.tet_bound.tri, ctx.hovered_facet);
			// FacetAttribute<int> flag(ctx.tet_bound.tri, -1);
			// algo::naive_tag({f}, flag, naive_constraints);
			// um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.tet_bound.tri, flag.ptr, "flag", ctx.mesh_);

			GEO::Attribute<GEO::signed_index_t> flag(
				ctx.mesh_.facets.attributes(), "flag"
			);

			algo::naive_tag({f}, flag, naive_constraints);

		}
		else if (current_algo == Smudge) {
			GEO::Attribute<GEO::signed_index_t> flag(
				ctx.mesh_.facets.attributes(), "flag"
			);

			if (smudge_ref_value == -1) {
				smudge_ref_value = flag[ctx.hovered_facet];
			}

			flag[ctx.hovered_facet] = smudge_ref_value;
		}
	}
}

void PaintFlagTool::paint_bucket() {
	// std::cout << "facet: " << ctx.selected_facet << " in group: " << facet_by_features[ctx.selected_facet] << std::endl;

	int feature = facet_by_features[ctx.selected_facet];

	if (current_algo == None) {
		GEO::Attribute<GEO::signed_index_t> flag(
			ctx.mesh_.facets.attributes(), "flag"
		);

		for (int i = 0; i < facet_by_features.size(); i++) {
			if (facet_by_features[i] != feature)
				continue;
			
			flag[i] = value;
		}
	}
	else if (current_algo == Naive) {
		// Extract facets indexes of selected feature
		std::vector<UM::Surface::Facet> facets;
		for (int i = 0; i < facet_by_features.size(); i++) {
			if (facet_by_features[i] != feature)
				continue;
			
			facets.push_back(UM::Surface::Facet(ctx.tet_bound.tri, i));
		}

		GEO::Attribute<GEO::signed_index_t> flag(
			ctx.mesh_.facets.attributes(), "flag"
		);

		algo::naive_tag(facets, flag, naive_constraints);

		// FacetAttribute<int> flag(ctx.tet_bound.tri, -1);
		// algo::naive_tag(facets, flag, naive_constraints);
		// um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(ctx.tet_bound.tri, flag.ptr, "flag", ctx.mesh_);
	}
}

void PaintFlagTool::attribute_transfert() {
	// Transfert attribute from surface tri to volume tet
	// Because it is volume tet attr that is used for polycubify (not the surface attr) !

	// Transfert "tri_flag" from UM to GEO
	FacetAttribute<int> tri_flag(ctx.tet_bound.tri, -1);
	CellFacetAttribute<int> tet_flag(ctx.tet_bound.tet, -1);
	um_bindings::um_attr_from_geo_attr<GEO::MESH_FACETS>(ctx.mesh_, "flag", ctx.tet_bound.tri, tri_flag.ptr);
	// Transfert from surface to volume
	ctx.tet_bound.set_attribute_to_volume(tri_flag, tet_flag);
	// Transfert "tet_flag" from GEO to UM
	um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(ctx.tet_bound.tet, tet_flag.ptr, "tet_flag", ctx.mesh_);
}

void PaintFlagTool::mouse_button_callback(int button, int action, int mods, int source) {

	if (ctx.is_facet_hovered() && current_mode == Charts) {
		// Paint
		paint_bucket();
	}
	
	// Transfert attr from surface to volume
	attribute_transfert();

	smudge_ref_value = -1;
}

void PaintFlagTool::scroll_callback(double xoffset, double yoffset) {}

void PaintFlagTool::validate_callback() {

}

bool PaintFlagTool::is_compatible() { 
	return !ctx.mesh_metadata.filename.empty() && ctx.mesh_metadata.cell_type == MESH_TET;
}

void PaintFlagTool::escape_callback() {

}

