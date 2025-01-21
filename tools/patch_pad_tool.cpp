#include "patch_pad_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"


bool PatchPadTool::draw_object_properties() {
	
	if (ImGui::Button("Compute features")) {

		// Clear existing feature lines
		ctx.mesh_.edges.clear();

		// Compute patches on hex boundary
		{
			DisjointSet ds(ctx.hex.nfacets());

			for (auto h : ctx.hex.iter_halfedges()) {
				
				if (!h.facet().on_boundary() || !h.opposite_f().opposite_c().active())
					continue;
				
				auto opp_f = h.opposite_f().opposite_c().opposite_f().facet();
				Quad3 q = h.facet();
				Quad3 opp_q = opp_f;

				// Compute angle
				if (abs(q.normal() * opp_q.normal()) > 0.6) {
					ds.merge(h.facet(), opp_f);
				}
			}
	
			ds.get_sets_id(patches);
			is_init_patches = true;
		}

		std::vector<std::pair<int, int>> edges;

		// Get patches edge borders
		for (auto h : ctx.hex.iter_halfedges()) {
			
			if (!h.facet().on_boundary())
				continue;

			auto opp_c = h.opposite_f().opposite_c();
			if (!opp_c.active()) {
				edges.push_back({h.from(), h.to()});
				// Add to edge border
				continue;
			}

			auto f = h.facet();
			auto opp_f = h.opposite_f().opposite_c().opposite_f().facet();
			
			if (patches[f] != patches[opp_f]) {
				// Add edge to border
				edges.push_back({h.from(), h.to()});
			}

		}

		// Add feature to model
		ctx.mesh_.edges.create_edges(edges.size());
		for (int e = 0; e < edges.size(); e++) {
			ctx.mesh_.edges.set_vertex(e, 0, edges[e].first);
			ctx.mesh_.edges.set_vertex(e, 1, edges[e].second);
		}

		// Compute boundary
		// HexBoundary hb(ctx.hex);

	}

	return true;
}

void PatchPadTool::draw_viewer_properties() {}

void PatchPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
		
}

void PatchPadTool::hover_callback(double x, double y, int source) {

}

void PatchPadTool::mouse_button_callback(int button, int action, int mods, int source) {
	
}

void PatchPadTool::scroll_callback(double xoffset, double yoffset) {

}

void PatchPadTool::validate_callback() {

}

bool PatchPadTool::is_compatible() { 
	return true;
}

