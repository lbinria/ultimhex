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

	if (ctx.gui_mode == BlocPadding || ctx.switch_mode == BlocPadding) {
		if (ImGui::Button("Facet selection")) {
			select_mode = 0;
		}
		if (ImGui::Button("Chart selection")) {
			compute_patches();
			select_mode = 1;
		}
	}


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
	if (/*step > 0 ||*/ !ctx.is_facet_hovered())
		return;


	std::vector<int> hovered_facets;

	// Chart selection mode
	if (select_mode == 1) {
		for (auto f : ctx.hex_bound->quad.iter_facets()) {
			if (patches[ctx.hex_bound->hex_facet(f)] == patches[ctx.hex_bound->hex_facet(ctx.hovered_facet)])
				hovered_facets.push_back(f);
		}
	} else {
		hovered_facets.push_back(ctx.hovered_facet);
	}

	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	// Reset last hovered facets as not hovered, if not selected
	for (auto last_hovered_f : last_hovered_facets) {
		if (last_hovered_f >= 0 && hovered_attr[last_hovered_f] == 1) {
			hovered_attr[last_hovered_f] = 0;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
		}
	}

	// Last hovered facet become current hovered facet
	last_hovered_facets = hovered_facets;

	for (auto f : hovered_facets) {
		// If not selected, facet is hovered
		if (hovered_attr[f] < 2) {
			hovered_attr[f] = 1;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 1;
		}

		if (ctx.left_mouse_pressed) {
			hovered_attr[f] = 2;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 2;
		} else if (ctx.right_mouse_pressed) {
			hovered_attr[f] = 0;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 0;
		}
	}

}

// void BlocPadTool::hover_callback(double x, double y, int source) {
// 	if (/*step > 0 ||*/ !ctx.is_facet_hovered())
// 		return;


// 	std::vector<int> hovered_facets;

// 	// Chart selection mode
// 	if (select_mode == 1) {
// 		for (auto f : ctx.hex_bound->quad.iter_facets()) {
// 			if (patches[f] == patches[ctx.hovered_facet])
// 				hovered_facets.push_back(f);
// 		}
// 	} else {
// 		hovered_facets.push_back(ctx.hovered_facet);
// 	}

// 	// Attribute hovered / selected, enable visualizing hovered / selected facets
// 	GEO::Attribute<int> hovered_attr(
// 		ctx.mesh_.facets.attributes(), "hovered"
// 	);
// 	GEO::Attribute<int> cell_facets_hovered_attr(
// 		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
// 	);

// 	// Reset last hovered facet as not hovered, if not selected
// 	if (last_hovered_f >= 0 && hovered_attr[last_hovered_f] == 1) {
// 		hovered_attr[last_hovered_f] = 0;
// 		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
// 	}

// 	// Last hovered facet become current hovered facet
// 	last_hovered_f = ctx.hovered_facet;

// 	// If not selected, facet is hovered
// 	if (hovered_attr[last_hovered_f] < 2) {
// 		hovered_attr[last_hovered_f] = 1;
// 		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 1;
// 	}

// 	if (ctx.left_mouse_pressed) {
// 		hovered_attr[last_hovered_f] = 2;
// 		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 2;
// 	} else if (ctx.right_mouse_pressed) {
// 		hovered_attr[last_hovered_f] = 0;
// 		cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
// 	}

// }

void BlocPadTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void BlocPadTool::scroll_callback(double xoffset, double yoffset) {

}

void BlocPadTool::compute_patches() {
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

void BlocPadTool::compute_feature_lines() {

	compute_patches();

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
}

void BlocPadTool::validate_callback() {

	// Attribute hovered / selected, enable visualizing hovered / selected facets
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
				
				if (!next_f.active())
					// Extremity is reached !
					cell_facets_hovered_attr[cur_f.opposite_f().next().next().opposite_f().facet()] = 2;

				cur_f = next_f;
			}


		}


		/*
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
		*/



		ctx.hex_bound = std::make_unique<HexBoundary>(ctx.hex, selected_cell);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

		// Attribute hovered / selected, enable visualizing hovered / selected facets
		GEO::Attribute<int> hovered_attr(
			ctx.mesh_.facets.attributes(), "hovered"
		);
		GEO::Attribute<int> cell_facets_hovered_attr(
			ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
		);

		for (auto f : ctx.hex_bound->quad.iter_facets()) {
			hovered_attr[f] = 2;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 2;
		}	



		ctx.view.change_mode(ViewBinding::Mode::Surface);
		

		



		// Compute feature lines
		// compute_feature_lines();


		step = 1;

	} else if (step == 1) {

		// Try to check validity
		bool valid = true;
		for (auto h : ctx.hex_bound->hex.iter_halfedges()) {
			if (cell_facets_hovered_attr[h.facet()] != 2)
				continue;


			// Check opp
			auto opp = h.opposite_f().opposite_c();
			if (!opp.active()) {
				continue;
			}

			int count = 0;
			for (auto eh : h.iter_CCW_around_edge()) {
				if (cell_facets_hovered_attr[eh.facet()] == 2)
					count++;
			}
			// if (cell_facets_hovered_attr[h.opposite_f().facet()] == 2)
			// 	count++;

			// opp = h.opposite_c();

			// if (opp.active() && cell_facets_hovered_attr[opp.opposite_f().facet()] == 2)
			// 	count++;

			// opp = h.opposite_f().opposite_c();

			// if (opp.active() && cell_facets_hovered_attr[opp.opposite_f().facet()] == 2)
			// 	count++;

			if (count != 2) {
				valid = false;
				break;
			}
		}

		std::cout << "Padding valid ? " << valid << std::endl;


		CellFacetAttribute<bool> pad_face(ctx.hex_bound->hex, false);
		
		GEO::Attribute<int> cell_facets_hovered_attr(
			ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
		);

		//
		for (auto f : ctx.hex_bound->hex.iter_facets()) {
			pad_face[f] = cell_facets_hovered_attr[f] == 2;
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

		// CellFacetAttribute<bool> pad_face(ctx.hex_bound->hex, false);
		
		// GEO::Attribute<bool> cell_filter(
		// 	ctx.mesh_.cells.attributes(), "filter"
		// );

		// //
		// for (auto f : ctx.hex_bound->hex.iter_facets()) {
		// 	if (!cell_filter[f.cell()])
		// 		continue;

		// 	pad_face[f] = !f.opposite().active() || !cell_filter[f.opposite().cell()];

		// 	// Exclude extremities ?
		// 	if (is_outgoing_padding && cell_facets_hovered_attr[f] == 2)
		// 		pad_face[f] = false;
		// }

		// // Beurk ! Need to encapsulate hex Bound !!!!
		// ctx.hex_bound = NULL;
		// BenjaminAPI::pad(ctx.hex, pad_face);

		// // Reconstruct
		// ctx.hex_bound = std::make_unique<HexBoundary>(ctx.hex);
		// um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);

		// ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		// switch_view();
		// ctx.view.change_mode(ViewBinding::Mode::Volume);
		// ctx.view.show_vertices_ = false;


	}


}

bool BlocPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void BlocPadTool::escape_callback() {
	// Reset mesh
	ctx.hex_bound = std::make_unique<HexBoundary>(ctx.hex_bound->hex);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

	// Clear tool
	clear();

	// Reset wireframe
	ctx.mesh_.edges.clear();

	// Return to default view
	switch_view();
}

