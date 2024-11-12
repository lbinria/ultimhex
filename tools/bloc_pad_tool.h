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


	bool is_compatible() override;
	void key_callback(int key, int scancode, int action, int mods) {}
	void escape_callback() override;

	void clear() override {
		step = 0;
		start_f_idx = -1;

		depth = 1;
		hovered_bloc_facets.clear();
		selected_bloc_facets.clear();
		hovered_bloc_cells.clear();
		selected_bloc_cells.clear();
		wireframe.clear();

		bloc_surface.clear();

		ctx.view.show_volume_ = true;
	}

	void switch_view(bool to_wireframe);

	int step = 0;
	int start_f_idx = -1;
	std::vector<int> selected_region_cells;

	int depth = 1;

	std::vector<int> hovered_bloc_facets;
	std::vector<int> selected_bloc_facets;
	std::vector<int> hovered_bloc_cells;
	std::vector<int> selected_bloc_cells;

	std::vector<int> bloc_surface;

	std::vector<UM::Segment3> wireframe;

};