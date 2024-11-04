#include "layer_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "geom_ultimaille_binding.h"
#include "context.h"
#include "gl_draw.h"

void loop_cut(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Facet&)> f) {
	assert(hex.connected());

	std::vector<bool> visited(24 * hex.ncells(), false);
	std::queue<int> q;
	q.push(start_he);
	visited[start_he] = true;


	while (!q.empty()) {

		auto he_idx = q.front();
		auto he = Volume::Halfedge(hex, he_idx);

		q.pop();

		auto opp_c = he.next().opposite_f().opposite_c();

		if (opp_c.active()) {
			auto n_he = opp_c.opposite_f().next();

			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		} else {
			auto n_he = he.next().opposite_f().next();
			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		}

		opp_c = he.prev().opposite_f().opposite_c();

		if (opp_c.active()) {
			auto n_he = opp_c.opposite_f().prev();

			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		} else {
			auto n_he = he.prev().opposite_f().prev();
			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		}

		// Propagate to front
		opp_c = he.opposite_f().next().next().opposite_f().opposite_c();

		if (opp_c.active()) {
			auto n_he = opp_c;

			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		} 

		// Process current
		auto facet = he.opposite_f().facet();
		f(facet);
		
	}
}

void loop_cut2(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge& /* cur_he */, bool/* on_border */)> f) {
	assert(hex.connected());
	std::vector<bool> visited(24 * hex.ncells(), false);

	auto cur_he = start_he;

	while (true) {

		f(cur_he, true);

		auto opp_c = cur_he.next().opposite_f().opposite_c();
		if (!opp_c.active() || cur_he.next().opposite_f().next().facet().on_boundary()) {
			cur_he = cur_he.next().opposite_f().next();

		} else {
			cur_he = opp_c.opposite_f().next();

			// on border ?
		}

		if (visited[cur_he])
			break;

		visited[cur_he] = true;


	}

}

void LayerPadTool::draw_gui() {
	if(ImGui::Button("Loop padding")) {
		ctx.gui_mode = LayerPadding;
		// TODO notify Change tool
		hovered_path.clear();
		selected_path.clear();
	}
}

void LayerPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color) {
	gl_draw::draw_path(hovered_path, hovered_color, true);
	gl_draw::draw_path(selected_path, selected_color, true);
}

void LayerPadTool::hover_callback(double x, double y, int source) {
	hovered_path.clear();

	if (ctx.hex.connected() && ctx.is_cell_hovered()) {

		Volume::Cell um_c(ctx.hex, ctx.hovered_cell);
		
		Volume::Halfedge hovered_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.hovered_edge, ctx.hovered_lfacet)));
		
		if (hovered_he.active()) {
			loop_cut2(ctx.hex, hovered_he, [&](Volume::Halfedge &he, bool on_border) {
				UM::vec3 a = he.from().pos();
				UM::vec3 b = he.to().pos();
				hovered_path.push_back(a);
				hovered_path.push_back(b);
			});
		}
	}
}

void LayerPadTool::mouse_button_callback(int button, int action, int mods, int source) {
	selected_path = hovered_path;
}

void LayerPadTool::validate_callback() {
	CellFacetAttribute<bool> pad_face(ctx.hex);

	Volume::Cell um_c(ctx.hex, ctx.selected_cell);
	Volume::Halfedge start_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.selected_edge, ctx.selected_lfacet)));

	loop_cut(ctx.hex, start_he, [&](UM::Volume::Facet &f) {
		pad_face[f] = true;
	});

	BenjaminAPI::pad(ctx.hex, pad_face);

	um_bindings::geo_mesh_from_um_hex(ctx.hex, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

	selected_path.clear();
	write_by_extension("padded.geogram", ctx.hex);
}

bool LayerPadTool::is_compatible() { 
	return ctx.mesh_metadata.cell_type == MESH_HEX;
}