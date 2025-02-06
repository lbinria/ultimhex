#include "patch_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"


bool PatchPadTool::draw_object_properties() {
	
	ImGui::Checkbox("Extends to concave", &extends_to_concave);

	if (ImGui::Button("Patch selection")) {
		compute_features();
		return true;
	}

	return true;
}

void PatchPadTool::compute_features() {
	// Clear existing feature lines
	ctx.mesh_.edges.clear();

	compute_patches_for_selection();

	// Compute feature lines
	// std::vector<std::pair<int, int>> edges;

	// // Get patches edge borders
	// int n1 = 0;
	// int n2 = 0;
	// int eee = 0;
	// for (auto h : ctx.hex_bound->quad.iter_halfedges()) {
	// 	eee++;
	// 	auto opp = h.opposite();

	// 	if (!opp.active()) {
	// 		// Add to edge border
	// 		edges.push_back({h.from(), h.to()});
	// 		n1++;
	// 		continue;
	// 	}

	// 	auto f = h.facet();
	// 	auto opp_f = opp.facet();
		
	// 	if (patches[f] != patches[opp_f]) {
	// 		// Add edge to border
	// 		edges.push_back({h.from(), h.to()});
	// 		n2++;
	// 	}

	// }

	// // Add feature to GEO model (allow to view them)
	// ctx.mesh_.edges.create_edges(edges.size());
	// for (int e = 0; e < edges.size(); e++) {
	// 	ctx.mesh_.edges.set_vertex(e, 0, edges[e].first);
	// 	ctx.mesh_.edges.set_vertex(e, 1, edges[e].second);
	// }

	// Display surface with feature lines
	ctx.view.change_mode(ViewBinding::Mode::Surface);
	ctx.view.attribute_ = "facets.hovered";
	ctx.view.attribute_name_ = "hovered";
	ctx.view.attribute_min_ = 0;
	ctx.view.attribute_max_ = 2;
	ctx.gui_mode = PatchPadding;
}

void PatchPadTool::compute_patches_for_selection2() {

	std::vector<int> border_h;
	
	DisjointSet ds(ctx.hex_bound->hex.nfacets());

	for (auto h : ctx.hex_bound->quad.iter_halfedges()) {
		
		if (!h.opposite().active()) {
			continue;
		}

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
		
		// We expect to merge adjacents facets that are almost coplanar
		if (angle > 0.4) {
			// ds.merge(f, opp_f);
			ds.merge(ctx.hex_bound->hex_facet(f), ctx.hex_bound->hex_facet(opp_f));
		} else {
			border_h.push_back(h);
		}

	}

	// Extends halfedge on border until that reach model border
	for (auto hi : border_h) {
		auto h = ctx.hex_bound->hex_halfedge(hi);
		
		auto cur_h = h;
		
		while (cur_h.opposite_f().opposite_c().active()) {
			auto opp = cur_h.opposite_f().opposite_c();
			ds.merge(cur_h.facet(), opp.opposite_f().facet());
			cur_h = opp.opposite_f().next().next();
		}
	}

	ds.get_sets_id(patches);
	is_init_patches = true;

}

void PatchPadTool::compute_patches_for_selection() {

	// border_h.clear();
	// Compute patches
	{
		DisjointSet ds(ctx.hex_bound->quad.nfacets());

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

			// We expect to merge adjacents facets that are almost coplanar
			bool should_merge = angle > 0.4;
			// We can also merge adjacents facets that are not coplanar but have a concave angle
			if (extends_to_concave)
				should_merge = should_merge ||  v > 0;

			if (should_merge) {
				ds.merge(f, opp_f);
			}

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


	int patch = patches[ctx.hovered_facet];

	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);

	{
		std::vector<int> old_hovered_facets = hovered_facets;

		// Remove old hovered facets as hovered
		for (auto f : old_hovered_facets) {
			hovered_attr[f] = 0;
		}
	}

	// Clear current hovered facet to recompute them
	hovered_facets.clear();
	
	// If model facet is the same patch as hovered facet, mark it as hovered
	for (auto f : ctx.hex_bound->quad.iter_facets()) {
		if (patches[f] != patch)
			continue;

		// Convert UM facet index to GEO facet index
		// auto gf = um_bindings::geo_facet_index_from_um_facet_index(f, 6);

		hovered_attr[f] = 1;
		hovered_facets.push_back(f);
	}

	// Keep selected facets as selected
	for (auto f : selected_facets) {
		hovered_attr[f] = 2;
	}
}

void PatchPadTool::mouse_button_callback(int button, int action, int mods, int source) {

	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);

	{
		std::vector<int> old_selected_facets = selected_facets;

		for (auto f : old_selected_facets) {
			hovered_attr[f] = 0;
		}

		for (auto f : hovered_facets) {
			hovered_attr[f] = 0;
		}
	}

	selected_facets = hovered_facets;

	for (auto f : selected_facets) {
		hovered_attr[f] = 2;
	}

}

void PatchPadTool::scroll_callback(double xoffset, double yoffset) {

}

void PatchPadTool::validate_callback() {

	CellFacetAttribute<bool> to_pad(ctx.hex_bound->hex, false);

	CellFacetAttribute<int> v_hovered(ctx.hex_bound->hex, -1);

	{
		FacetAttribute<int> s_hovered(ctx.hex_bound->quad, -1);
		um_bindings::um_attr_from_geo_attr<GEO::MESH_FACETS>(ctx.mesh_, "hovered", ctx.hex_bound->quad, s_hovered.ptr);
		// Transfert attribute from surface to volume
		ctx.hex_bound->set_attribute_to_volume(s_hovered, v_hovered);
	}

	for (auto f : ctx.hex_bound->hex.iter_facets()) {
		// If facet is selected, pad it
		if (v_hovered[f] == 2) {
			to_pad[f] = true;
		}
	}

	// // Extends halfedge on border until that reach model border
	// for (auto hi : border_h) {
	// 	auto h = ctx.hex_bound->hex_halfedge(hi);
		
	// 	auto cur_h = h;
		
	// 	while (cur_h.opposite_f().opposite_c().active()) {
	// 		auto opp = cur_h.opposite_f().opposite_c();
	// 		cur_h = opp.opposite_f().next().next();
	// 		to_pad[cur_h.facet()] = true;
	// 	}
	// }

	// Necessary for altering hex, because hex / quad share points !
	ctx.hex_bound->clear_surface();
	BenjaminAPI::pad(ctx.hex_bound->hex, to_pad);


	// Clear tool
	clear();
	// Clear computed patches (as they should be recomputed)
	clear_patches();

	// Recreate surface from hex
	// Should update hex bound surface
	ctx.hex_bound = std::make_unique<HexBoundary>(ctx.hex);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

	ctx.view.change_mode(ViewBinding::Mode::Volume);

	// write_by_extension("patch_pad.geogram", ctx.hex_bound->hex, {});
	// write_by_extension("patch_pad_surf.geogram", ctx.hex_bound->quad, {});
}

bool PatchPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

