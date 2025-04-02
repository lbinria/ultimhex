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
#include <nicostuff/algo/framework/benjamin_API.h>

#include "geom_ultimaille_binding.h"

#include "mesh_metadata.h"

// std libs
#include <optional>

enum GUIMode {
	Camera = 0,
	Hover = 1,
	Filter = 2,
	Painting = 3,
	LayerPadding = 4,
	BlocPadding = 5,
	PatchPadding = 6,
	NewBlocPadding = 7,
	Polycubify = 8,
	HexCollapse = 9,
	Smooth = 10,
	Embedit = 11,
	PathConstraintPadding = 12
};

struct ViewBinding {

	bool &show_attributes_;
	bool &show_vertices_;
	bool &show_volume_;
	bool &show_surface_;

	bool &show_hexes_;

	enum Mode {
		None,
		Surface,
		Volume
	};

	const char * modes[3] = { "None", "Surface", "Volume" };
	int current_mode = Surface;


	index_t &current_colormap_index_;
	std::vector<unsigned int> colormaps;
	std::string &attribute_;
	GEO::MeshElementsFlags &attribute_subelements_;
	std::string &attribute_name_;
	float &attribute_min_;
	float &attribute_max_;
	bool &lighting_;

	void change_mode(Mode mode) {
		current_mode = mode;
		if (mode == ViewBinding::Mode::Surface) {
			show_surface_ = true;
			show_volume_ = false;
		} else if (mode == ViewBinding::Mode::Volume) {
			show_surface_ = false;
			show_volume_ = true;
		} else {
			show_surface_ = false;
			show_volume_ = false;
		}
	}

	void change_mode(int i) {
		change_mode((ViewBinding::Mode)i);
	}

	void switch_to_volume_select_mode() {
		change_mode(ViewBinding::Mode::Volume);
		attribute_subelements_ = GEO::MeshElementsFlags::MESH_CELLS;
		attribute_ = "cells.cell_hovered";
		attribute_name_ = "cell_hovered";
		attribute_min_ = 0;
		attribute_max_ = 2;
	}

	void switch_to_surface_select_mode() {
		change_mode(ViewBinding::Mode::Surface);
		attribute_subelements_ = GEO::MeshElementsFlags::MESH_FACETS;
		attribute_ = "facets.hovered";
		attribute_name_ = "hovered";
		attribute_min_ = 0;
		attribute_max_ = 2;
	}

	ImVec4 getColorFromTexture1D(GLuint textureID, int x) {
		
		glBindTexture(GL_TEXTURE_2D, textureID);

		// Get the width
		GLint width;
		glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width);

		unsigned char pixels[width * 4];
		glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
		
		int i = x * 4;	

		return ImVec4(pixels[i] / 255.f, pixels[i + 1] / 255.f, pixels[i + 2] / 255.f, pixels[i + 3]);
	}


};

struct Context {

	Context(Mesh &mesh, MeshGfx &mesh_gfx, ViewBinding view) : mesh_(mesh), mesh_gfx_(mesh_gfx), view(view) {
		hex_bound = std::make_unique<MyHexBoundary>(hex);
		tet_bound = std::make_unique<TetBoundary>(tet);
	}

	bool show_last_picked_point_ = false;
	bool show_hovered_cell_overlay_ = true;
	bool show_hovered_cell_facet_overlay_ = false;
	float overlay_thickness = 3.;

	bool region_selection_activated = false;

	index_t hovered_vertex = NO_VERTEX;
	index_t hovered_edge = NO_EDGE;
	index_t hovered_facet = NO_FACET;
	index_t hovered_cell_facet = NO_FACET;
	index_t hovered_cell_lfacet = NO_FACET;
	index_t hovered_cell = NO_CELL;

	index_t last_hovered_vertex = NO_VERTEX;
	index_t last_hovered_edge = NO_EDGE;
	index_t last_hovered_facet = NO_FACET;
	index_t last_hovered_cell_facet = NO_FACET;
	index_t last_hovered_cell_lfacet = NO_FACET;
	index_t last_hovered_cell = NO_CELL;


	std::set<index_t> hovered_cells;

	index_t selected_vertex = NO_VERTEX;
	index_t selected_edge = NO_EDGE;
	index_t selected_facet = NO_FACET;
	index_t selected_cell_facet = NO_FACET;
	index_t selected_cell_lfacet = NO_FACET;
	index_t selected_cell = NO_CELL;
	std::vector<index_t> selected_cells;

	int brush_size = 1;

	int um_hovered_cell_facet() { 
		int n_facet_per_cell =mesh_metadata.cell_type == MESH_HEX ? 6 : 4;
		return um_bindings::um_facet_index_from_geo_facet_index(hovered_cell_facet, n_facet_per_cell);
	}

	int hovered_lhe() {
		return um_bindings::he_from_cell_e_lf(hovered_edge, hovered_cell_lfacet);
	}

	Volume::Halfedge hovered_he(Volume &v) {
		Volume::Cell c(v, hovered_cell);
		return c.halfedge(hovered_lhe());
	}

	// TODO Please refactor this !
	bool is_cell_hovered() { return (hovered_cell != NO_CELL && hovered_cell < mesh_.cells.nb()); }
	bool is_facet_hovered() { return (hovered_facet != NO_FACET && hovered_facet < mesh_.facets.nb()); }
	bool is_cell_facet_hovered() { return (hovered_cell_facet != NO_FACET && hovered_cell_facet < mesh_.cell_facets.nb()); }
	bool is_cell_lfacet_hovered() { return (hovered_cell_lfacet != NO_FACET && hovered_cell_lfacet < 6); }
	bool is_cell_edge_hovered() { return (hovered_edge != NO_EDGE && hovered_edge < 12); }

	bool is_cell_selected() { return (selected_cell != NO_CELL && selected_cell < mesh_.cells.nb()); }
	bool is_facet_selected() { return (selected_facet != NO_FACET && selected_facet < mesh_.facets.nb()); }
	bool is_cell_facet_selected() { return (selected_cell_facet != NO_FACET && selected_cell_facet < mesh_.cell_facets.nb()); }
	bool is_cell_lfacet_selected() { return (selected_cell_lfacet != NO_FACET && selected_cell_lfacet < 6); }
	bool is_cell_edge_selected() { return (selected_edge != NO_EDGE && selected_edge < 12); }

	void reset_hovered_selected() {
		hovered_vertex = NO_VERTEX;
		hovered_edge = NO_EDGE;
		hovered_cell_facet = NO_FACET;
		hovered_cell_lfacet = NO_FACET;
		hovered_cell = NO_CELL;
		selected_vertex = NO_VERTEX;
		selected_edge = NO_EDGE;
		selected_cell_facet = NO_FACET;
		selected_cell_lfacet = NO_FACET;
		selected_cell = NO_CELL;

		hovered_cells.clear();
		selected_cells.clear();
		brush_size = 1;
	}


	void recompute_hex() {
		// Recompute hex boundary
		hex_bound = std::make_unique<MyHexBoundary>(hex);
		
		// Make chart segmentation & init embedding
		tri_chart = std::make_unique<FacetAttribute<int>>(tet_bound->tri, -1);
		quad_chart = std::make_unique<FacetAttribute<int>>(hex_bound->quad, -1);
		BenjaminAPI::embeditinit(tet_bound->tri, *tri_chart, hex_bound->hex, *emb_attr, *quad_chart, false);

		// Recompute geo mesh
		um_bindings::geo_mesh_from_hexboundary(*hex_bound, mesh_);
		// Refresh view
		mesh_gfx_.set_mesh(&mesh_);
	}



	GEO::vec3 click_pos;
	bool left_mouse_pressed = false;
	bool right_mouse_pressed = false;
	GUIMode gui_mode = Camera;
	GUIMode switch_mode = Camera;

	GEO::vec3 posA;
	GEO::vec3 posB;

	Mesh &mesh_;
	MeshGfx &mesh_gfx_;

	UM::Tetrahedra tet;
	UM::Hexahedra hex;

	std::unique_ptr<TetBoundary> tet_bound;
	std::unique_ptr<MyHexBoundary> hex_bound;

	// Embedding
	std::unique_ptr<CellFacetAttribute<int>> emb_attr;
	std::unique_ptr<FacetAttribute<int>> tri_chart;
	std::unique_ptr<FacetAttribute<int>> quad_chart;
	std::vector<int> emb;

	MeshMetadata mesh_metadata;



	ViewBinding view;
};