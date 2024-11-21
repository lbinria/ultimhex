#pragma once

#include <ultimaille/all.h>
#include <geogram_gfx/third_party/imgui/imgui.h>

#include "tool.h"


struct PaintFlagTool : public Tool {

	PaintFlagTool(Context &ctx) : Tool(ctx) {}

	std::string get_name() { return "Flag painting"; }

	void init();

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
		value = -1;
		current_mode = Facet;
		current_algo = None;
		smudge_ref_value = -1;
		facet_by_features.clear();
		facet_by_color.clear();
		for (int i = 0; i < 3; i++)
			naive_constraints[i] = true;
	}

	void paint_bucket();
	void attribute_transfert();

	enum PaintMode {
		Facet,
		Charts
	};
	enum PaintAlgo {
		None,
		Smudge,
		Naive
	};

	const char * modes[2] = { "Paint facets", "Paint charts" };
	const char * algos[3] = { "None", "Smudge", "Naive" };
	int current_mode = Facet;
	int current_algo = None;
	int value = -1;

	int smudge_ref_value = -1;

	bool naive_constraints[3] = {true, true, true};

	std::vector<int> facet_by_features;
	std::vector<int> facet_by_color;



};