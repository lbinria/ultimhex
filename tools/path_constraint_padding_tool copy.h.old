#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct PathConstraintPaddingTool : public Tool {

	PathConstraintPaddingTool(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Path constraint padding"; }

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
		layer.clear();
		last_hovered_facets.clear();
		v_cell_facets_hovered_attr.clear();
		mode = Select;
	}

	void reset_attr();
	void puff_view(Volume::Halfedge &h);
	void slice();

	std::vector<int> layer;
	std::vector<int> last_hovered_facets;

	GEO::vector<int> v_cell_facets_hovered_attr;

	float v_plane_center[3] = {0., 0., 0.};
	
	float v_plane_rot_x = 0.;
	float v_plane_rot_y = 0.;
	float v_plane_rot_z = 0.;
	GEO::vec3 plane_center;
	GEO::vec3 plane_rot;
	GEO::vec3 v_plane_corners[4];

	enum Mode {
		Select,
		Preview
	};

	Mode mode = Select;

};