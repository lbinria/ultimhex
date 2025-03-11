#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>
#include <nicostuff/algo/framework/benjamin_API.h>

#include "tool.h"


struct HexCollapseTool : public Tool {

	HexCollapseTool(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Hex collapse"; }

	bool draw_object_properties() override;
	void draw_viewer_properties() override;

	void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) override; 
	void mouse_button_callback(int button, int action, int mods, int source) override;
	void scroll_callback(double xoffset, double yoffset) override;
	void hover_callback(double x, double y, int source) override;
	void validate_callback() override;


	void compute_patches();
	void compute_feature_lines();

	bool is_compatible() override;
	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() override;


	void compute_layers();

	void clear() override {
		layers.clear();
		last_hovered_cells.clear();
		hovered_cells.clear();
		selected_cells.clear();
		selected_layer = -1;
		hovered_h = -1;
	}

	void reset() {
		clear();
		compute_layers();
	}

	private:

	std::vector<int> layers;
	std::vector<int> last_hovered_cells;
	std::vector<int> hovered_cells;
	std::vector<int> selected_cells;
	int selected_layer = -1;
	int hovered_h = -1;


};