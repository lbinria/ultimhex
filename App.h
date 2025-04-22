#pragma once

#include <geogram_gfx/gui/application.h>                // for set_style()
#include <geogram_gfx/gui/simple_mesh_application.h>
#include <geogram/basic/file_system.h>                  // for is_file(), extension() in load()
#include <geogram/basic/command_line.h>                 // for get_arg_bool() in load()
#include <geogram/mesh/mesh_io.h>                       // for MeshIOFlags, mesh_load() in load()
#include <geogram/basic/command_line.h>                 // for CmdLine::get_arg() and CmdLine::set_arg()
#include <geogram/basic/string.h>                       // for String::string_ends_with()
#include "gui_base.h"

#include <ultimaille/all.h>

#include "geom_ultimaille_binding.h"

#include "context.h"
#include "mesh_metadata.h"

// Tools
#include "tools/filter_tool.h"
#include "tools/camera_tool.h"
#include "tools/hover_tool.h"
#include "tools/layer_stack_redefinition_tool.h"
#include "tools/paint_flag_tool.h"
#include "tools/bloc_pad_tool.h"
#include "tools/patch_pad_tool.h"
#include "tools/polycubify_tool.h"
#include "tools/hex_collapse_tool.h"
#include "tools/smooth_tool.h"
#include "tools/embedit_tool.h"
#include "tools/path_constraint_padding_tool.h"
#include "tools/hex_split_tool.h"

// std libs
#include <optional>

// new colormaps
#define COLORMAP_HOVER_SELECTION (SIMPLE_APPLICATION_NB_COLORMAPS)
#define COLORMAP_FLAGGING (SIMPLE_APPLICATION_NB_COLORMAPS + 1)


class App : public SimpleMeshApplicationExt {

public:

    App(const std::string name = "ultimhex");

	void ImGui_initialize() override;

protected:

	void normalize_mesh();
	void view_mesh();

	bool refresh_hovered();

	void draw_scene() override;

	void draw_gui() override;

	void draw_menu_bar() override;

	void draw_viewer_properties() override;

    void draw_object_properties() override;

    bool load(const std::string& filename) override;

	std::string supported_write_file_extensions() override;
	std::string supported_read_file_extensions() override;

	bool save(const std::string& filename) override;

	void GL_initialize() override;

	void mouse_button_callback(int button, int action, int mods, int source) override;

	void scroll_callback(double xoffset, double yoffset) override;

	void cursor_pos_callback(double x, double y, int source) override;

	void key_callback(int key, int scancode, int action, int mods) override;

	void labeling_visu_mode_transition();

	void reset();

protected:

	ColorArray hover_selection_colors_; // colors for hover & selection
	ColorArray flagging_colors_; // colors for flagging

	// Declare tools
	CameraTool camera_tool;
	HoverTool hover_tool;
	FilterTool filter_tool;
	PaintFlagTool paint_flag_tool;
	BlocPadTool bloc_pad_tool;
	PatchPadTool patch_pad_tool;
	LayerPad2 new_bloc_pad_tool;
	PolycubifyTool polycubify_tool;
	HexCollapseTool hex_collapse_tool;
	SmoothTool smooth_tool;
	EmbeditTool embedit_tool;
	PathConstraintPaddingTool path_constraint_padding_tool;
	HexSplitTool hex_split_tool;


	// std::size_t nb_tools = std::size(GUIMode::Camera);
	std::unique_ptr<Tool> tools[13] = {
		std::make_unique<CameraTool>(camera_tool), 
		std::make_unique<HoverTool>(hover_tool), 
		std::make_unique<FilterTool>(filter_tool), 
		std::make_unique<PaintFlagTool>(paint_flag_tool), 
		std::make_unique<BlocPadTool>(bloc_pad_tool),
		std::make_unique<PatchPadTool>(patch_pad_tool),
		std::make_unique<LayerPad2>(new_bloc_pad_tool),
		std::make_unique<PolycubifyTool>(polycubify_tool),
		std::make_unique<HexCollapseTool>(hex_collapse_tool),
		std::make_unique<SmoothTool>(smooth_tool),
		std::make_unique<EmbeditTool>(embedit_tool),
		std::make_unique<PathConstraintPaddingTool>(path_constraint_padding_tool),
		std::make_unique<HexSplitTool>(hex_split_tool),
	};

	enum MeshVolumeType {
		Tet,
		Hex
	};

	bool is_loading = false;

	bool load_step = false;
	// float size_factor = 0.2;
	float size_factor = 0.5;
	char * gmsh_path = "gmsh";

	// TODO move to gui_base
	bool show_grid_ = true;
	bool show_axes_ = true;
	bool show_features = true;

	int shrink = 0;


	// Option
	bool tool_preview = true;

	// List of current meshes
	std::map<std::string, Mesh> meshes;



	// TODO rename
	Context context_;

	UM::vec2 last_mouse_pos;

};