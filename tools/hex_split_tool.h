#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>
#include <nicostuff/algo/framework/benjamin_API.h>


#include "../helpers.h"
#include "tool.h"
#include "layer_selector.h"

struct HexSplitTool : public Tool {

	HexSplitTool(Context &ctx) : Tool(ctx), layer_selector(ctx) {}

	std::string get_name() { return "Hex split"; }

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



	void clear() override {
		layer_selector.clear();
	}

	void reset();

	private:
	LayerSelector layer_selector;


	bool auto_smooth = false;


};