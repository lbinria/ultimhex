#include "bloc_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"


std::vector<int> extract_region_3d(UM::Hexahedra &hex, std::vector<int> facets, int depth) {

	std::vector<int> cells;

	for (int i = 0; i < depth; i++) {

		std::vector<int> next_facets;

		for (auto f_idx : facets) {
			Volume::Facet f(hex, f_idx);

			auto opp_c = f.halfedge(0).opposite_f().next().next().opposite_f().opposite_c();
			if (opp_c.active()) {
				next_facets.push_back(opp_c.facet());
			}

			cells.push_back(f.cell());
		}

		facets = next_facets;

		if (next_facets.empty())
			break;
	}

	return cells;
}

// /**
//  * Search a rect region between two facets
//  */
// std::optional<std::vector<int>> extract_region_between_facets(UM::Hexahedra &hex, UM::Volume::Volume::Facet &start_f, UM::Volume::Facet &end_f) {
// 	assert(hex.connected());

// 	std::vector<bool> from_to_visited(24 * hex.ncells(), false);
// 	std::vector<bool> to_from_visited(24 * hex.ncells(), false);

// 	std::queue<int> q;

// 	int dir = 0;
	
// 	q.push(start_f);
// 	from_to_visited[start_f] = true;

// 	while (!q.empty()) {

// 		auto f_idx = q.front();
// 		Volume::Facet f(hex, f_idx);
	
// 		q.pop();

// 		// Found end facet
// 		if (f == end_f) {
// 			break;
// 		}

// 		for (int ij = 0; ij < 2; ij++) {

// 			auto h = f.halfedge((dir + ij) % 4);

// 			auto opp_c = h.next().opposite_f().opposite_c();
// 			if (opp_c.active()) {
// 				auto n_h = opp_c.opposite_f().next();
// 				auto n_f = n_h.facet();

// 				if (!from_to_visited[n_f]) {
// 					q.push(n_f);
// 					from_to_visited[n_f] = true;
// 				}
// 			}

// 		}

// 	}




// 	std::vector<int> facets;

// 	for (int i = 0; i < from_to_visited.size(); i++) {
// 		if (!from_to_visited[i])
// 			continue;

// 		Volume::Facet f(hex, i);
// 		facets.push_back(f);
// 	}
// 	return facets;

// 	// return std::nullopt;
// }

/**
 * Search a rect region between two facets
 */
std::optional<std::vector<int>> extract_region_between_facets(UM::Hexahedra &hex, UM::Volume::Volume::Facet &start_f, UM::Volume::Facet &end_f) {
	assert(hex.connected());

	// std::vector<bool> visited(24 * hex.ncells(), false);

	// Choose a direction {+1,+1}, {-1,+1}, {-1,-1}, {+1,-1}
	for (int dir = 0; dir < 4; dir++) {

		std::vector<std::tuple<int, int, int>> facets_by_coords;

		auto x_h = start_f.halfedge(dir);
		auto y_h = start_f.halfedge((dir + 1) % 4);

		auto cur_x_h = x_h;
		auto cur_y_h = y_h;

		
		int dy = 0;
		while (true) {

			int dx = 0;
			while (true) {
				
				// Found !
				auto facet = cur_x_h.facet();
				facets_by_coords.push_back(std::tuple(facet, dx, dy));

				if (facet == end_f) {
					// Should remove facets that are out of rect
					std::vector<int> facets;
					for (auto facet_by_coords : facets_by_coords) {
						if (std::get<1>(facet_by_coords) <= dx && std::get<2>(facet_by_coords) <= dy) {
							facets.push_back(std::get<0>(facet_by_coords));
						}
					}
					// std::cout << "coords: " << dx << ", " << dy << std::endl;
					return facets;
				}

				auto opp_c = cur_x_h.opposite_f().opposite_c();
				// if (!opp_c.active() || visited[cur_x_h])
				if (!opp_c.active()) {
					break;
				}

				// visited[cur_x_h] = true;
				auto nxt_x_he = opp_c.opposite_f().next().next();

				cur_x_h = nxt_x_he;
				dx++;

			}

			auto opp_c = cur_y_h.opposite_f().opposite_c();
			// if (!opp_c.active() || visited[cur_y_h]) {
			if (!opp_c.active()) {
				break;
			}

			// visited[cur_y_h] = true;
			cur_y_h = opp_c.opposite_f().next().next();
			cur_x_h = cur_y_h.facet().halfedge(dir);

			dy++;
		}

	}

	return std::nullopt;
}

/**
 * Extract facets from the surface of a bloc of cells
 */
std::vector<int> extract_surf_facet(Hexahedra &hex, std::vector<int> &cells) {
	assert(hex.connected());

	std::vector<int> facets;

	for (auto c_idx : cells) {
		UM::Volume::Cell c(hex, c_idx);
		for (auto f : c.iter_facets()) {
			// Contains
			bool is_in_set = f.opposite().active() ? std::find(cells.begin(), cells.end(), f.opposite().cell()) != cells.end() : false;

			// if (f.on_boundary())
			if (f.on_boundary() || !is_in_set)
				facets.push_back(f);
		}
	}

	return facets;
}

std::vector<UM::Segment3> extract_wireframe(Hexahedra &hex) {

	std::vector<UM::Segment3> segments;

	HexBoundary hex_bound(hex);

	for (auto h : hex_bound.quad.iter_halfedges()) {

		Volume::Halfedge hh(hex, hex_bound.hex_halfedge(h));
		if (hh.opposite_f().opposite_c().active())
			continue;

		
		Segment3 s = h;
		segments.push_back(s);
	}

	return segments;
}

bool BlocPadTool::draw_object_properties() {
	if(ImGui::Button("Bloc selection")) {
		ctx.gui_mode = BlocPadding;
		return true;
	}

	if (step == 2) {
		if (ImGui::InputInt("depth", &depth, 1, 1)) {
			selected_bloc_cells = extract_region_3d(ctx.hex, selected_bloc_facets, depth);
		}
	}

	return false;
}

void BlocPadTool::draw_viewer_properties() {}

void BlocPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	
	// glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(1.0f, 1.0f, 0.3f, 1.0f).data());
	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(1.f,0.f,0.f).data());

	glupBegin(GLUP_LINES);

	for (auto s : wireframe) {
		glupPrivateVertex3dv(um_bindings::geo_vec(s.a).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(s.b).data());
	}

	glupEnd();

	for (auto c : hovered_bloc_cells) {
		gl_draw::draw_cell_overlay(ctx.mesh_, c, colorMapInfo, 0.0, ctx.overlay_thickness);
	}
	
	for (auto c : selected_bloc_cells) {
		gl_draw::draw_cell_overlay(ctx.mesh_, c, colorMapInfo, 0.5, ctx.overlay_thickness);
	}

	// for (auto f_idx : bloc_surface) {
	// 	auto f = Volume::Facet(ctx.hex, f_idx);
	// 	auto c = f.cell();
	// 	auto lf = f.id_in_cell();
	// 	gl_draw::draw_cell_facet_overlay(ctx.mesh_, c, lf, colorMapInfo, 0.5, 5);
	// }



	
}

void BlocPadTool::hover_callback(double x, double y, int source) {

	if (!(ctx.hex.connected() && ctx.is_cell_hovered() && ctx.is_cell_facet_hovered()))
		return;

	if (step == 0) {
	
		hovered_bloc_cells.clear();
		hovered_bloc_cells.push_back(ctx.hovered_cell);
		hovered_bloc_facets.clear();
		hovered_bloc_facets.push_back(ctx.hovered_cell_facet);

	} else if (step == 1) {
		// Get UM facets TODO compute this directly in hover to avoid call umbinding each time ! 
		Volume::Facet start_f(ctx.hex, start_f_idx);
		Volume::Facet end_f(ctx.hex, um_bindings::um_facet_index_from_geo_facet_index(ctx.hovered_cell_facet, 6));

		// Get selection from start facet to end facet
		auto facets_opt = extract_region_between_facets(ctx.hex, start_f, end_f);

		if (facets_opt.has_value()) {
			// std::cout << "coords:" << coord_opt.value().first << ", " << coord_opt.value().second << std::endl;
			hovered_bloc_cells.clear();
			hovered_bloc_facets.clear();

			for (auto f : facets_opt.value()) {
				int c = Volume::Facet(ctx.hex, f).cell();
				hovered_bloc_cells.push_back(c);
				hovered_bloc_facets.push_back(f);
			}
		}
	}
}

void BlocPadTool::mouse_button_callback(int button, int action, int mods, int source) {
	
	if (ctx.hex.connected() && ctx.is_cell_hovered() && ctx.is_cell_facet_hovered()) {

		if (step == 0 || step >= 2) {
			if (step >= 2) {
				clear();
			}

			start_f_idx = Volume::Facet(ctx.hex, um_bindings::um_facet_index_from_geo_facet_index(ctx.selected_cell_facet, 6));
			step++;

		} else if (step == 1) {
			selected_bloc_cells = hovered_bloc_cells;
			selected_bloc_facets = hovered_bloc_facets;
			step++;
		}

	}
}

void BlocPadTool::scroll_callback(double xoffset, double yoffset) {
	if (step < 2)
		return;

	if (yoffset > 0) {
		depth++;
	} else if (yoffset < 0 && depth > 1) {
		depth--;
	}

	// // Exit wireframe view
	// if (depth == 1 && wireframe.size() > 0) {
	// 	switch_view(false);
	// }

	if (yoffset != 0 && depth >= 1) {
		// Enter wireframe view
		if (wireframe.size() == 0)
			switch_view(true);
		
		selected_bloc_cells = extract_region_3d(ctx.hex, selected_bloc_facets, depth);
		// bloc_surface = extract_surf_facet(ctx.hex, selected_bloc_cells);

	}

}

void BlocPadTool::switch_view(bool to_wireframe) {
	if (to_wireframe) {
		// Extract wireframe and hide volume + surface
		wireframe = extract_wireframe(ctx.hex);
		ctx.view.change_mode(ViewBinding::Mode::None);
	} else {
		wireframe.clear();
		ctx.view.change_mode(ViewBinding::Mode::Volume);
	}
}

void BlocPadTool::validate_callback() {

	// Can only validate at step 2
	if (step != 2)
		return;

	// Extract facets on the surface of selection
	auto facets = extract_surf_facet(ctx.hex, selected_bloc_cells);
	// tmp_facets = facets;

	// Set as facet to pad
	CellFacetAttribute<bool> pad_face(ctx.hex);
	for (auto f : facets) {
		pad_face[f] = true;
	}

	// Pad !
	BenjaminAPI::pad(ctx.hex, pad_face);

	// Refresh GEO mesh from UM mesh
	um_bindings::geo_mesh_from_um_hex(ctx.hex, ctx.mesh_);



	// Refresh view	
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	write_by_extension("padded.geogram", ctx.hex);

	// Clear tool
	clear();



	// TEST !

	// ctx.mesh_gfx_.set_cells_color(MESH_HEX, 1.0, 1.0, 1.0, 0.5);
	// HexBoundary hex_bound(ctx.hex);
	// um_bindings::geo_mesh_from_hexboundary(hex_bound, ctx.mesh_);
	// ctx.mesh_gfx_.set_surface_color(1.0,1.0,1.0,0.2);

	// GEO::Attribute<bool> filter(
	// 	ctx.mesh_.cells.attributes(), "filter"
	// );
	// for (auto f : facets) {
	// 	auto c = Volume::Facet(ctx.hex, f).cell();
	// 	filter[c] = true;
	// }
	// ctx.mesh_gfx_.set_filter(MESH_CELLS, "filter");


	// wireframe = extract_wireframe(ctx.hex);
	// std::cout << wireframe.size() << std::endl;
	// std::cout << ctx.hex.ncorners() << std::endl;
}

bool BlocPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void BlocPadTool::escape_callback() {
	// Clear tool
	clear();
}

