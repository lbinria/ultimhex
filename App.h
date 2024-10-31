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

#include "mesh_metadata.h"

// std libs
#include <optional>

// new colormaps
#define COLORMAP_HOVER_SELECTION (SIMPLE_APPLICATION_NB_COLORMAPS)


class App : public SimpleMeshApplicationExt {
public:

    App(const std::string name = "polycubify");

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


	enum MeshVolumeType {
		Tet,
		Hex
	};

	enum GUIMode {
		Camera,
		Hover,
		Painting,
		LoopPadding,
		BlocPadding
	};

	bool is_loading = false;

	// TODO move to gui_base
	bool show_grid_ = true;
	bool show_axes_ = true;

	bool show_last_picked_point_ = false;
	bool show_hovered_cell_overlay_ = true;
	bool show_hovered_cell_facet_overlay_ = false;
	float overlay_thickness = 3.;

	// Option
	bool tool_preview = true;

	// List of current meshes
	std::map<std::string, Mesh> meshes;

	// TODO rename all to hovered_sometin
	index_t hovered_vertex = NO_VERTEX;
	index_t hovered_edge = NO_EDGE;
	// TODO rename to cellfacet
	index_t hovered_facet = NO_FACET;
	index_t hovered_lfacet = NO_FACET;
	index_t hovered_cell = NO_CELL;

	index_t selected_vertex = NO_VERTEX;
	index_t selected_edge = NO_EDGE;
	index_t selected_facet = NO_FACET;
	index_t selected_lfacet = NO_FACET;
	index_t selected_cell = NO_CELL;

	// TODO Please refactor this !
	bool is_cell_hovered() { return (hovered_cell != NO_CELL && hovered_cell < mesh_.cells.nb()); }
	bool is_cell_facet_hovered() { return (hovered_facet != NO_FACET && hovered_facet < mesh_.cell_facets.nb()); }
	bool is_cell_lfacet_hovered() { return (hovered_lfacet != NO_FACET && hovered_lfacet < 6); }
	bool is_cell_edge_hovered() { return (hovered_edge != NO_EDGE && hovered_edge < 12); }

	bool is_cell_selected() { return (selected_cell != NO_CELL && selected_cell < mesh_.cells.nb()); }
	bool is_cell_facet_selected() { return (selected_facet != NO_FACET && selected_facet < mesh_.cell_facets.nb()); }
	bool is_cell_lfacet_selected() { return (selected_lfacet != NO_FACET && selected_lfacet < 6); }
	bool is_cell_edge_selected() { return (selected_edge != NO_EDGE && selected_edge < 12); }

	void reset_hovered_selected() {
		hovered_vertex = NO_VERTEX;
		hovered_edge = NO_EDGE;
		hovered_facet = NO_FACET;
		hovered_lfacet = NO_FACET;
		hovered_cell = NO_CELL;
		selected_vertex = NO_VERTEX;
		selected_edge = NO_EDGE;
		selected_facet = NO_FACET;
		selected_lfacet = NO_FACET;
		selected_cell = NO_CELL;
	}

	GEO::vec3 click_pos;
	bool left_mouse_pressed = false;
	GUIMode gui_mode = Camera;
	GUIMode switch_mode = Camera;
	int paint_value = 0;

	GEO::vec3 posA;
	GEO::vec3 posB;

	// Path on surface mesh to display
	std::vector<UM::vec3> hovered_path;
	std::vector<UM::vec3> selected_path;



	int bloc_pad_step = 0;
	// int bloc_pad_0_selection = -1;
	// int bloc_pad_1_selection = -1;
	// std::vector<cell_cross> blocpad0[4];
	// std::vector<cell_cross> blocpad1[4];
	int bloc_start_f = -1;
	int bloc_end_f = -1;
	std::vector<int> hovered_bloc_cells;
	std::vector<int> selected_bloc_cells;


	// UM::vec3 posAb;
	// UM::vec3 posBb;
	// UM::vec3 posN;


	UM::Tetrahedra tet;
	UM::Hexahedra hex;
	TetBoundary tet_bound;

	MeshMetadata mesh_metadata;

	std::vector<std::pair<int, UM::vec3>> flag_dirs;

	// int he_n = 0;
	// um_bindings::MeshBinding mesh_binding;

	std::vector<UM::vec3> um_facetus;

};