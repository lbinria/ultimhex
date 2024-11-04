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


// std libs
#include <optional>

#include "layer_pad_tool.h"
#include "bloc_pad_tool.h"

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

	// Bloc padding data (move to tool)
	struct cell_cross {
		index_t c;
		index_t lf;
		int dist;
	};

	// std::optional<std::tuple<int, int, int>> App::compute_bloc_pad_intersection(std::vector<cell_cross> cc0s[4], std::vector<cell_cross> cc1s[4]);


protected:

	ColorArray hover_selection_colors_; // colors for hover & selection


	LayerPadTool layer_pad_tool;
	BlocPadTool bloc_pad_tool;

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



	Context context_;


	// UM::vec3 posAb;
	// UM::vec3 posBb;
	// UM::vec3 posN;




	std::vector<std::pair<int, UM::vec3>> flag_dirs;

	// int he_n = 0;
	// um_bindings::MeshBinding mesh_binding;

	std::vector<UM::vec3> um_facetus;

};