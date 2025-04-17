#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct EmbeditTool : public Tool {

	EmbeditTool(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Embedit"; }

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
		n_charts = 0;
		selected_chart_color = 0;
		mode = Paint;
		view_mode = Quad;
		is_init = false;
	}

	enum Mode {
		Pick,
		Paint
	};

	enum ViewMode {
		Tri,
		Quad
	};

	Mode mode;
	ViewMode view_mode = Quad;

	bool auto_smooth = false;
	bool is_init = false;
	int n_charts = 0;
	int selected_chart_color = 0;


};