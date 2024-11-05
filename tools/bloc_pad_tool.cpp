#include "bloc_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

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
					std::cout << "coords: " << dx << ", " << dy << std::endl;
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

std::vector<int> extract_surf_facet(Hexahedra &hex, std::vector<int> &cells) {
	assert(hex.connected());

	std::vector<int> facets;
	for (auto c_idx : cells) {
		UM::Volume::Cell c(hex, c_idx);
		for (auto f : c.iter_facets()) {
			if (!f.on_boundary() || f.opposite().active())
				continue;
			
			facets.push_back(f);
		}
	}

	return facets;
}

bool BlocPadTool::draw_gui() {
	if(ImGui::Button("Bloc padding")) {
		ctx.gui_mode = BlocPadding;
		return true;
	}

	return false;
}

void BlocPadTool::draw_viewer_properties() {}

void BlocPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	for (auto c : hovered_bloc_cells) {
		gl_draw::draw_cell_overlay(ctx.mesh_, c, colorMapInfo, 0.0, ctx.overlay_thickness);
	}
	for (auto c : selected_bloc_cells) {
		gl_draw::draw_cell_overlay(ctx.mesh_, c, colorMapInfo, 0.5, ctx.overlay_thickness);
	}
}

void BlocPadTool::hover_callback(double x, double y, int source) {

	if (bloc_pad_step != 1)
		return;

	if (ctx.hex.connected() && ctx.is_cell_facet_hovered()) {

		Volume::Facet um_f(ctx.hex, um_bindings::um_facet_index_from_geo_facet_index(ctx.hovered_facet, 6));
		Volume::Facet bloc_start_ff(ctx.hex, bloc_start_f);
		auto facets_opt = extract_region_between_facets(ctx.hex, bloc_start_ff, um_f);
		if (facets_opt.has_value()) {
			// std::cout << "coords:" << coord_opt.value().first << ", " << coord_opt.value().second << std::endl;
			hovered_bloc_cells.clear();
			for (auto f : facets_opt.value()) {
				int c = Volume::Facet(ctx.hex, f).cell();
				hovered_bloc_cells.push_back(c);
			}
		}
	}
}

void BlocPadTool::mouse_button_callback(int button, int action, int mods, int source) {
	
	if (ctx.hex.connected() && ctx.is_cell_facet_hovered()) {

		if (bloc_pad_step == 0) {
			Volume::Facet um_f(ctx.hex, um_bindings::um_facet_index_from_geo_facet_index(ctx.selected_facet, 6));
			bloc_start_f = um_f;
			bloc_pad_step++;
		}
		else if (bloc_pad_step == 1) {
			selected_bloc_cells = hovered_bloc_cells;
			bloc_pad_step++;
		}

	}
}

void BlocPadTool::validate_callback() {

	if (bloc_pad_step != 2)
		return;


	auto facets = extract_surf_facet(ctx.hex, selected_bloc_cells);
	CellFacetAttribute<bool> pad_face(ctx.hex);
	for (auto f : facets) {
		pad_face[f] = true;
	}
	BenjaminAPI::pad(ctx.hex, pad_face);

	um_bindings::geo_mesh_from_um_hex(ctx.hex, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

	
	write_by_extension("padded.geogram", ctx.hex);

	// Clear tool
	selected_bloc_cells.clear();
	hovered_bloc_cells.clear();
	bloc_pad_step = 0;
	bloc_start_f = -1;
}

bool BlocPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}

void BlocPadTool::escape_callback() {
	// Clear tool
	selected_bloc_cells.clear();
	hovered_bloc_cells.clear();
	bloc_pad_step = 0;
	bloc_start_f = -1;
}

