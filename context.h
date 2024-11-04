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

enum GUIMode {
	Camera,
	Hover,
	Painting,
	LayerPadding,
	BlocPadding
};

struct Context {

	Context(Mesh &mesh, MeshGfx &mesh_gfx) : mesh_(mesh), mesh_gfx_(mesh_gfx), tet_bound(tet) {}

	bool show_last_picked_point_ = false;
	bool show_hovered_cell_overlay_ = true;
	bool show_hovered_cell_facet_overlay_ = false;
	float overlay_thickness = 3.;

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

	Mesh &mesh_;
	MeshGfx &mesh_gfx_;

	UM::Tetrahedra tet;
	UM::Hexahedra hex;
	TetBoundary tet_bound;

	MeshMetadata mesh_metadata;
};