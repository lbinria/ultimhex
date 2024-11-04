#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct LayerPadTool : public Tool {

	LayerPadTool(Context &ctx) : Tool(ctx) {}

	static constexpr char* mode = "LayerPadding";

	bool draw_gui() override;
	void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) override; 
	void mouse_button_callback(int button, int action, int mods, int source) override;
	void hover_callback(double x, double y, int source) override;
	void validate_callback() override;


	bool is_compatible() override;
	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() {}

	void clear() override {
		// Clear path
		hovered_path.clear();
		selected_path.clear();
	}

	// Path on surface mesh to display
	std::vector<UM::vec3> hovered_path;
	std::vector<UM::vec3> selected_path;

};