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
#include <camera_tool.h>
#include <hover_tool.h>
#include <paint_flag_tool.h>
#include "tools/layer_pad_tool.h"
#include "tools/bloc_pad_tool.h"
#include <polycubify_tool.h>

// std libs
#include <optional>

// new colormaps
#define COLORMAP_HOVER_SELECTION (SIMPLE_APPLICATION_NB_COLORMAPS)


class App : public SimpleMeshApplicationExt {

public:

    App(const std::string name = "ultimhex");

	void ImGui_initialize() override;

protected:

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

	void cursor_pos_callback(double x, double y, int source) override;

	void key_callback(int key, int scancode, int action, int mods) override;

	void labeling_visu_mode_transition();

	void reset();

protected:

	ColorArray hover_selection_colors_; // colors for hover & selection

	// Declare tools
	CameraTool camera_tool;
	HoverTool hover_tool;
	PaintFlagTool paint_flag_tool;
	LayerPadTool layer_pad_tool;
	BlocPadTool bloc_pad_tool;
	PolycubifyTool polycubify_tool;

	// std::size_t nb_tools = std::size(GUIMode::Camera);
	std::unique_ptr<Tool> tools[6] = {
		std::make_unique<CameraTool>(camera_tool), 
		std::make_unique<HoverTool>(hover_tool), 
		std::make_unique<PaintFlagTool>(paint_flag_tool), 
		std::make_unique<LayerPadTool>(layer_pad_tool), 
		std::make_unique<BlocPadTool>(bloc_pad_tool),
		std::make_unique<PolycubifyTool>(polycubify_tool)
	};

	enum MeshVolumeType {
		Tet,
		Hex
	};

	bool is_loading = false;

	// TODO move to gui_base
	bool show_grid_ = true;
	bool show_axes_ = true;



	// Option
	bool tool_preview = true;

	// List of current meshes
	std::map<std::string, Mesh> meshes;



	// TODO rename
	Context context_;

	// TOOD remove
	// UM::vec3 posAb;
	// UM::vec3 posBb;
	// UM::vec3 posN;
	// int he_n = 0;

	std::vector<std::pair<int, UM::vec3>> flag_dirs;

};