#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"

struct PatchPadTool : public Tool {

	PatchPadTool(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Patch pad tools"; }

	bool draw_object_properties() override;
	void draw_viewer_properties() override;
	void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) override; 
	void mouse_button_callback(int button, int action, int mods, int source) override;
	void scroll_callback(double xoffset, double yoffset) override;
	void hover_callback(double x, double y, int source) override;
	void validate_callback() override;


	void fix_padding();

	bool is_compatible() override;

	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() {
		clear();
	}

	void clear() override {
		
		// Attribute hovered / selected, enable visualizing hovered / selected facets
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

	void clear_patches() {
		patches.clear();
		is_init_patches = false;
	}

	void compute_patches_for_selection();
	void compute_patches_for_selection2();
	
	void extend_patches();
	void compute_features();

	float threshold = 80.f;
	bool extends_to_concave = false;
	bool traverse = false;

	bool is_init_patches = false;
	std::vector<int> patches;


};