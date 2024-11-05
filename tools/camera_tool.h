#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct CameraTool : public Tool {

	CameraTool(Context &ctx) : Tool(ctx) {}

	static constexpr char* mode = "Camera";

	bool draw_gui() override;
	void draw_viewer_properties() override;

	void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) override; 
	void mouse_button_callback(int button, int action, int mods, int source) override;
	void hover_callback(double x, double y, int source) override;
	void validate_callback() override;


	bool is_compatible() override;
	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() override;

	void clear() override {
		bloc_pad_step = 0;
		bloc_start_f = -1;
		bloc_end_f = -1;
		hovered_bloc_cells.clear();
		selected_bloc_cells.clear();
	}

	int bloc_pad_step = 0;
	int bloc_start_f = -1;
	int bloc_end_f = -1;
	std::vector<int> hovered_bloc_cells;
	std::vector<int> selected_bloc_cells;

};