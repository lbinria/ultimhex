#include "patch_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"


bool PatchPadTool::draw_object_properties() {
	
	ImGui::Checkbox("Extends to concave", &extends_to_concave);

	ImGui::SliderFloat("Max angle", &threshold, 0.f, 180.f);
	
	if (ImGui::Button("Puff pastry 0##btn_patch_pad_tool_puff_0", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		std::vector<int> layer;
		helpers::get_halfedge_layers(ctx.hex_bound->hex, layer);

		CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
		helpers::puff(ctx.hex_bound->hex, *ctx.hex_bound->hex.iter_halfedges().begin(), layer, selected);

		ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		// ctx.recompute_hex(selected);
		ctx.view.switch_to_surface_select_mode();
		is_puff_view = true;

	}

	if (ImGui::Button("Puff pastry 1##btn_patch_pad_tool_puff_1", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		std::vector<int> layer;
		helpers::get_halfedge_layers(ctx.hex_bound->hex, layer);

		CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
		helpers::puff(ctx.hex_bound->hex, (*ctx.hex_bound->hex.iter_halfedges().begin()).opposite_f().prev(), layer, selected);

		ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		// ctx.recompute_hex(selected);
		ctx.view.switch_to_surface_select_mode();
		is_puff_view = true;

	}

	if (ImGui::Button("Puff pastry 2##btn_patch_pad_tool_puff_2", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		std::vector<int> layer;
		helpers::get_halfedge_layers(ctx.hex_bound->hex, layer);

		CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
		helpers::puff(ctx.hex_bound->hex, (*ctx.hex_bound->hex.iter_halfedges().begin()).opposite_f().prev().opposite_f().prev(), layer, selected);

		ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		// ctx.recompute_hex(selected);
		ctx.view.switch_to_surface_select_mode();
		is_puff_view = true;

	}


	if (ImGui::Button("Patch selection")) {
		compute_features();
		// Display surface with feature lines
		ctx.view.switch_to_surface_select_mode();
		ctx.gui_mode = GUIMode::PatchPadding;
		return true;
	}

	return false;
}

void PatchPadTool::compute_features() {
	// Clear existing feature lines
	ctx.mesh_.edges.clear();
	compute_patches_for_selection();
}

void PatchPadTool::compute_patches_for_selection() {

	// Compute patches
	{
		DisjointSet ds(ctx.hex_bound->hex.nfacets());

		for (auto h : ctx.hex_bound->quad.iter_halfedges()) {
			
			if (!h.opposite().active())
				continue;

			auto f = h.facet();
			auto opp_f = h.opposite().facet();

			Quad3 q = f;
			Quad3 opp_q = opp_f;

			// Compute angle
			double angle = q.normal() * opp_q.normal();

			// Compute signed volume between 2 tetrahedrons of quad facets to found convex angles
			UM::vec3 bary;
			UM::vec3 bary_opp;
			for (auto qh : f.iter_halfedges()) {
				bary += qh.from().pos();
			}
			bary /= 4.;
			for (auto qh : opp_f.iter_halfedges()) {
				bary_opp += qh.from().pos();
			}
			bary_opp /= 4.;

			double v = Tetrahedron(h.from().pos(), h.to().pos(), bary, bary_opp).volume();


			double th = cos(threshold * (M_PI / 180.));

			// We expect to merge adjacents facets that are almost coplanar
			// if (angle > 0.4) {
			// if (angle > 0.95 || v > 0) {
			if (angle > th) {
				auto hf = ctx.hex_bound->hex_facet(f);
				auto opp_hf = ctx.hex_bound->hex_facet(opp_f);
				ds.merge(hf, opp_hf);
			} /* else if (v > 0) {

				// Extends halfedge on border until that reach model border
				auto cur_h = ctx.hex_bound->hex_halfedge(h);
				
				while (cur_h.opposite_f().opposite_c().active()) {
					auto opp = cur_h.opposite_f().opposite_c();
					ds.merge(cur_h.facet(), opp.opposite_f().facet());
					cur_h = opp.opposite_f().next().next();
				}

			}*/




		}

		ds.get_sets_id(patches);
		is_init_patches = true;
	}
}

void PatchPadTool::draw_viewer_properties() {}

void PatchPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
		
}

void PatchPadTool::hover_callback(double x, double y, int source) {

	if (!ctx.is_facet_hovered() || !is_init_patches || patches.empty())
		return;


	auto hex_hovered_f = ctx.hex_bound->hex_facet(ctx.hovered_facet);
	int patch = patches[hex_hovered_f];

	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	// Reset, clear old hovered facets, keep selected facets
	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		if (cell_facets_hovered_attr[f] == 2)
			continue;

		if (patches[f] == patch)
			cell_facets_hovered_attr[f] = 1;
		else
			cell_facets_hovered_attr[f] = 0;

		if (ctx.hex_bound->quad_facet(f) >= 0)
			hovered_attr[ctx.hex_bound->quad_facet(f)] = cell_facets_hovered_attr[f];
	}

}

void PatchPadTool::mouse_button_callback(int button, int action, int mods, int source) {

	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	// Remove all previous selected facets
	for (auto f : ctx.hex_bound->hex.iter_facets()) {

		if (cell_facets_hovered_attr[f] == 2)
			cell_facets_hovered_attr[f] = 0;

		if (ctx.hex_bound->quad_facet(f) >= 0)
			hovered_attr[ctx.hex_bound->quad_facet(f)] = cell_facets_hovered_attr[f];
	}

	// Hovered facets become selected facets
	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		if (cell_facets_hovered_attr[f] == 1)
			cell_facets_hovered_attr[f] = 2;
		
		if (ctx.hex_bound->quad_facet(f) >= 0)
			hovered_attr[ctx.hex_bound->quad_facet(f)] = cell_facets_hovered_attr[f];
	}

}

void PatchPadTool::scroll_callback(double xoffset, double yoffset) {

}

void PatchPadTool::validate_callback() {

	CellFacetAttribute<bool> to_pad(ctx.hex_bound->hex, false);

	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		if (cell_facets_hovered_attr[f] == 2) {
			to_pad[f] = true;
		}
	}

	// Necessary for altering hex, because hex / quad share points !
	ctx.hex_bound->clear_surface();
	BenjaminAPI::pad(ctx.hex_bound->hex, to_pad, *ctx.emb_attr);


	// Clear tool
	clear();
	// Clear computed patches (as they should be recomputed)
	clear_patches();

	// Recreate surface from hex
	// Should update hex bound surface
	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

	ctx.view.change_mode(ViewBinding::Mode::Volume);
	

	{
		// Visualize embedding
		FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
		ctx.hex_bound->set_attribute_to_surface(*ctx.emb_attr, surf_emb_attr);
		write_by_extension("after_patch_pad_emb.geogram", ctx.hex_bound->quad, {{}, {{"emb", surf_emb_attr.ptr}}, {}});
	}
}

bool PatchPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

