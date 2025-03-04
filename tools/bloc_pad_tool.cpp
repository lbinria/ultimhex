#include "bloc_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"


bool BlocPadTool::draw_object_properties() {

	if (ImGui::Button("Bloc selection")) {

		clear();
		switch_view();
		ctx.gui_mode = BlocPadding;

		return true;
	}

	ImGui::Checkbox("Outgoing padding ?", &is_outgoing_padding);

	return false;
}

void BlocPadTool::switch_view() {
	ctx.view.change_mode(ViewBinding::Mode::Surface);
	ctx.view.attribute_ = "facets.hovered";
	ctx.view.attribute_name_ = "hovered";
	ctx.view.attribute_min_ = 0;
	ctx.view.attribute_max_ = 2;
	ctx.mesh_gfx_.unset_filters();
}

void BlocPadTool::draw_viewer_properties() {}

void BlocPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
}

void BlocPadTool::hover_callback(double x, double y, int source) {
	if (step > 0 || !ctx.is_facet_hovered())
		return;


	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	// Set last hovered facet as not hovered, if not selected
	if (last_hovered_f >= 0 && hovered_attr[last_hovered_f] < 2) {
		hovered_attr[last_hovered_f] = 0;
		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
	}

	// Last hovered facet become current hovered facet
	last_hovered_f = ctx.hovered_facet;

	// If not selected, facet is hovered
	if (hovered_attr[last_hovered_f] < 2) {
		hovered_attr[last_hovered_f] = 1;
		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 1;
	}

	if (ctx.left_mouse_pressed) {
		hovered_attr[last_hovered_f] = 2;
		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 2;
	} else if (ctx.right_mouse_pressed) {
		hovered_attr[last_hovered_f] = 0;
		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
	}

}

void BlocPadTool::mouse_button_callback(int button, int action, int mods, int source) {

	if (step == 1 && ctx.is_cell_hovered()) {
		GEO::Attribute<bool> cell_filter(
			ctx.mesh_.cells.attributes(), "filter"
		);

		cell_filter[ctx.hovered_cell] = false;
		ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_CELLS);

	}

	// if (step > 0 && last_hovered_f < 0)
	// 	return;

}

void BlocPadTool::scroll_callback(double xoffset, double yoffset) {

}

void BlocPadTool::compute_patches_for_selection() {

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

			// Try to make traversing padding
			bool success = false;

			// double th = cos(threshold * (M_PI / 180.));

			// We expect to merge adjacents facets that are almost coplanar
			// if (angle > 0.4) {
			// if (angle > 0.95 || v > 0) {
			if (angle > 0.4) {
				auto hf = ctx.hex_bound->hex_facet(f);
				auto opp_hf = ctx.hex_bound->hex_facet(opp_f);
				ds.merge(hf, opp_hf);
			}

			ds.get_sets_id(patches);
			// is_init_patches = true;
		}
	}
}

void BlocPadTool::validate_callback() {

	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	if (step == 0) {

		int n_cells = 0;
		CellAttribute<bool> selected_cell(ctx.hex_bound->hex, false);

		// Extends facet selection to cells lace
		for (auto f : ctx.hex_bound->hex.iter_facets()) {

			if (cell_facets_hovered_attr[f] != 2)
				continue;

			// Extrude
			auto cur_f = f.halfedge(0);
			while (cur_f.active()) {
				
				if (!selected_cell[cur_f.cell()]) {
					selected_cell[cur_f.cell()] = true;
					n_cells++;
				}

				// Get opp facets
				auto next_f = cur_f.opposite_f().next().next().opposite_f().opposite_c();
				
				if (!next_f.active() || cell_facets_hovered_attr[next_f] == 3 /* reach a stop facet */)
					// Extremity is reached !
					cell_facets_hovered_attr[cur_f.opposite_f().next().next().opposite_f().facet()] = 2;

				cur_f = next_f;
			}


		}


		// Filter
		GEO::Attribute<bool> cell_filter(
			ctx.mesh_.cells.attributes(), "filter"
		);


		for (auto c : ctx.hex_bound->hex.iter_cells()) {
			cell_filter[c] = selected_cell[c];
		}

		// Switch view to preview selection
		ctx.view.change_mode(ViewBinding::Mode::Volume);
		ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_CELLS);



		// Compute feature lines
		compute_patches_for_selection();

		std::vector<std::pair<int, int>> edges;

		// Get patches edge borders
		int n1 = 0;
		int n2 = 0;
		int eee = 0;
		for (auto h : ctx.hex_bound->quad.iter_halfedges()) {
			eee++;
			auto opp = h.opposite();

			if (!opp.active()) {
				// Add to edge border
				edges.push_back({h.from(), h.to()});
				n1++;
				continue;
			}

			auto f = h.facet();
			auto opp_f = opp.facet();
			
			if (patches[f] != patches[opp_f]) {
				// Add edge to border
				edges.push_back({h.from(), h.to()});
				n2++;
			}

		}

		ctx.mesh_.edges.clear();
		// Add feature to GEO model (allow to view them)
		ctx.mesh_.edges.create_edges(edges.size());
		for (int e = 0; e < edges.size(); e++) {
			ctx.mesh_.edges.set_vertex(e, 0, edges[e].first);
			ctx.mesh_.edges.set_vertex(e, 1, edges[e].second);
		}

		step = 1;

	} else if (step == 1) {

		ctx.hex_preview.connect();

		CellFacetAttribute<bool> pad_face(ctx.hex_bound->hex, false);
		
		GEO::Attribute<bool> cell_filter(
			ctx.mesh_.cells.attributes(), "filter"
		);

		//
		for (auto f : ctx.hex_bound->hex.iter_facets()) {
			if (!cell_filter[f.cell()])
				continue;

			pad_face[f] = !f.opposite().active() || !cell_filter[f.opposite().cell()];

			// Exclude extremities ?
			if (is_outgoing_padding && cell_facets_hovered_attr[f] == 2)
				pad_face[f] = false;
		}

		// Beurk ! Need to encapsulate hex Bound !!!!
		ctx.hex_bound = NULL;
		BenjaminAPI::pad(ctx.hex, pad_face);

		// Reconstruct
		ctx.hex_bound = std::make_unique<HexBoundary>(ctx.hex);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);

		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		switch_view();
		ctx.view.change_mode(ViewBinding::Mode::Volume);
		ctx.view.show_vertices_ = false;


	}


}

bool BlocPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void BlocPadTool::escape_callback() {
	// Clear tool
	clear();

	// Reset wireframe
	ctx.mesh_.edges.clear();



	switch_view();

}

