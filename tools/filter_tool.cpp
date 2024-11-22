#include "filter_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"

bool FilterTool::draw_object_properties() {

	// if (ImGui::CollapsingHeader("Filter")) {
		
		GEO::Attribute<bool> filter(
			ctx.mesh_.facets.attributes(), "filter"
		);

		if (ImGui::Button("Reset filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {

			for (auto f : ctx.mesh_.facets)
				filter[f] = true;

			ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_FACETS);
		}

		if (ImGui::Button("Switch filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {

			for (auto f : ctx.mesh_.facets)
				filter[f] = !filter[f];

			ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_FACETS);
		}

		if (ImGui::Button("Chart filter", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
			ctx.gui_mode = Filter;
			mode = Chart;
		}
	// }

	return false;
}

void FilterTool::draw_viewer_properties() {

}

void FilterTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void FilterTool::hover_callback(double x, double y, int source) {

}

void FilterTool::mouse_button_callback(int button, int action, int mods, int source) {
	
	if (ctx.is_facet_hovered() && mode == Chart) {
		filter_chart();
	}

}

void FilterTool::filter_chart() {


	std::vector<int> facet_by_features;

	// TODO duplicate code in paint flag tool
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

	int feature = facet_by_features[ctx.selected_facet];

	GEO::Attribute<bool> filter(
		ctx.mesh_.facets.attributes(), "filter"
	);

	for (int i = 0; i < facet_by_features.size(); i++) {
		filter[i] = facet_by_features[i] != feature;
	}

	ctx.mesh_gfx_.set_filter(GEO::MeshElementsFlags::MESH_FACETS);
}

void FilterTool::scroll_callback(double xoffset, double yoffset) {}

void FilterTool::validate_callback() {

}

void FilterTool::escape_callback() {
	
}

bool FilterTool::is_compatible() { 
	return true;
}