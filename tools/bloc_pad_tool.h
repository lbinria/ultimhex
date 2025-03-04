#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct BlocPadTool : public Tool {

	BlocPadTool(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Bloc padding"; }

	bool draw_object_properties() override;
	void draw_viewer_properties() override;

	void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) override; 
	void mouse_button_callback(int button, int action, int mods, int source) override;
	void scroll_callback(double xoffset, double yoffset) override;
	void hover_callback(double x, double y, int source) override;
	void validate_callback() override;


	void compute_patches_for_selection();

	bool is_compatible() override;
	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() override;


	void switch_view();

	void clear() override {
		step = 0;
		last_hovered_f = -1;
		// hexs.clear();
		// hex_preview_facet_2_hex_facet.clear();
		patches.clear();

		// Clear selection attribute
		GEO::Attribute<int> hovered_attr(
			ctx.mesh_.facets.attributes(), "hovered"
		);
		GEO::Attribute<int> cell_facets_hovered_attr(
			ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
		);
		
		for (auto f : ctx.hex_bound->quad.iter_facets())
			hovered_attr[f] = 0;
		for (auto f : ctx.hex_bound->hex.iter_facets())
			cell_facets_hovered_attr[f] = 0;
	}

	int step = 0;
	int last_hovered_f;
	bool is_outgoing_padding = true;

	std::vector<int> patches;



};