#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct LayerPad2 : public Tool {

	LayerPad2(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Redefine layers stack"; }

	bool draw_object_properties() override;
	void draw_viewer_properties() override;

	void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) override; 
	void mouse_button_callback(int button, int action, int mods, int source) override;
	void scroll_callback(double xoffset, double yoffset) override;
	void hover_callback(double x, double y, int source) override;
	void validate_callback() override;


	bool is_compatible() override;
	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() override;

	void clear() override {
		is_init_layers = false;
		selected_h_idx = -1;

		hovered_cells.clear();
		selected_cells.clear();
		selected_cell_facets.clear();
		layers.clear();
	}

	// Pre-computation of layers in hex mesh
	bool is_init_layers = false;
	std::vector<int> layers;

	int selected_h_idx;
	std::vector<int> hovered_cells;
	std::vector<int> selected_cells;
	std::vector<std::pair<int, int>> selected_cell_facets;
	int hh;

	int n_layers_requested = 3;
	

};