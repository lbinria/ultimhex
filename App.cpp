// Internal libs
#include <ultimaille/all.h>
#include <geogram_gfx/full_screen_effects/ambient_occlusion.h>
#include <nicostuff/algo/framework/benjamin_API.h>

#include "App.h"

#include "geom_ultimaille_binding.h"
#include "tag_face.h"

#include "my_tet_boundary.h"
#include "gl_draw.h"
#include "mesh_metadata.h"

// JSON !!
#include <json.hpp>

// std libs
#include <queue>

#include "layer_pad_tool.h"

using json = nlohmann::json;

App::App(const std::string name) : 
	SimpleMeshApplicationExt(name), 
	context_(mesh_, mesh_gfx_),
	hover_selection_colors_({
		{1.f,0.3f,0.6f, 1.f}, // pink hover
		{0.95f,0.2f,0.0f, 1.f}  // red select
	}),
	layer_pad_tool(context_),
	bloc_pad_tool(context_)
{

}

void App::ImGui_initialize() {
    Application::ImGui_initialize();
    set_style("Dark");
    if(GEO::FileSystem::is_file("gui.ini")) {
        // Layout modification, saved with ImGui::SaveIniSettingsToDisk()
        // Larger docked object properties panel
        ImGui::LoadIniSettingsFromDisk("gui.ini");
    }
}

void App::draw_scene() {
    SimpleMeshApplicationExt::draw_scene();

	auto hovered_color = GEO::vec4f(1.f,0.3f,0.6f, 1.f);
	auto selected_color = GEO::vec4f(1.f,0.2f,0.0f, 1.f);

	// GLUPfloat *view;
	// glupGetMatrixfv(GLUP_MODELVIEW_MATRIX, view);
	// ImGuizmo::ViewManipulate(view, 2.f, ImVec2(100,100), ImVec2(20,20), ImU32(100));

	glupSetPointSize(10.0);

	if (show_axes_)
		gl_draw::draw_axis();
	if (show_grid_)
		gl_draw::draw_grid();

	auto hover_selection_colormap = colormaps_[COLORMAP_HOVER_SELECTION];

	// PATH
	if (context_.gui_mode == LayerPadding) {
		layer_pad_tool.draw(hovered_color, selected_color, hover_selection_colormap);
	} 
	if (context_.gui_mode == BlocPadding) {
		bloc_pad_tool.draw(hovered_color, selected_color, hover_selection_colormap);
	}

	// Overlays 
	// Cell facet
	if (context_.gui_mode == Hover) {
		if (context_.show_hovered_cell_facet_overlay_) {
			if (context_.is_cell_hovered() && context_.is_cell_facet_hovered()) {
				gl_draw::draw_cell_facet_overlay(mesh_, context_.hovered_cell, context_.hovered_lfacet, colormaps_[COLORMAP_HOVER_SELECTION], 0.0, context_.overlay_thickness);	
			}
		}
		// Cell
		if (context_.show_hovered_cell_overlay_) {
			if (context_.is_cell_hovered()) {
				gl_draw::draw_cell_overlay(mesh_, context_.hovered_cell, colormaps_[COLORMAP_HOVER_SELECTION], 0.0, context_.overlay_thickness);
			}
		}
		if (context_.is_cell_selected() && context_.is_cell_lfacet_selected()) {
			gl_draw::draw_cell_facet_overlay(mesh_, context_.selected_cell, context_.selected_lfacet, colormaps_[COLORMAP_HOVER_SELECTION], 0.5);
		}
		if (context_.is_cell_selected()) {
			gl_draw::draw_cell_overlay(mesh_, context_.selected_cell, colormaps_[COLORMAP_HOVER_SELECTION], 0.5);
		}
	}



	// Last picked point position as point
	if (context_.show_last_picked_point_) {
		glupBegin(GLUP_POINTS);
		glupPrivateVertex3dv(context_.click_pos.data());
		glupEnd();
	}

	for (auto p : um_facetus) {
		glupBegin(GLUP_POINTS);
		glupPrivateVertex3dv(um_bindings::geo_vec(p).data());
		glupEnd();
	}

	// Just test

	for (auto x : flag_dirs) {
		int flag = x.first;
		auto p = x.second;

		UM::vec3 dir{0,0,0};
		// UGLY but just for testing !
		if (flag == 0)
			dir = {-1,0,0};
		else if (flag == 1)
			dir = {0,-1,0};
		else if (flag == 2)
			dir = {0,0,-1};
		else if (flag == 3)
			dir = {1,0,0};
		else if (flag == 4)
			dir = {0,1,0};
		else if (flag == 5)
			dir = {0,0,1};

		gl_draw::draw_arrow(p, (p + dir * 0.05), 0.01, 8, 0.75, GEO::vec4f(1,0,0,1));
	}

}

void App::draw_gui() {
    SimpleMeshApplicationExt::draw_gui();
}

void App::draw_menu_bar() {
    SimpleApplication::draw_menu_bar();

    if(ImGui::BeginMainMenuBar()) {
        ImGui::EndMainMenuBar();
    }
}

void App::draw_viewer_properties() {
    SimpleMeshApplicationExt::draw_viewer_properties();

	ImGui::Checkbox("Show grid", &show_grid_);
	ImGui::Checkbox("Show axes", &show_axes_);
	ImGui::Separator();
	ImGui::Checkbox("Show vertices", &show_vertices_);
	ImGui::Checkbox("Show surface", &show_surface_);
	ImGui::Checkbox("Show volume", &show_volume_);
	ImGui::Separator();
	ImGui::Checkbox("Show picket point", &context_.show_last_picked_point_);
	ImGui::Separator();
	if (context_.gui_mode == Hover) {
		ImGui::Checkbox("Show cell overlay", &context_.show_hovered_cell_overlay_);
		ImGui::Checkbox("Show cell facet overlay", &context_.show_hovered_cell_facet_overlay_);
		ImGui::SliderFloat("Thickness", &context_.overlay_thickness, 1., 5.);
	}

}

void App::GL_initialize() {
    SimpleMeshApplicationExt::GL_initialize();
	init_rgba_colormap("hover_selection", 2, 1, hover_selection_colors_.as_chars());
    // state_transition(state_); // not all state_transition() code has been executed if GL was not initialized (in particular because missing colormaps)
}

std::vector<std::pair<int, UM::vec3>> compute_patches(UM::Triangles &tri, FacetAttribute<int> &tri_flag) {

	DisjointSet ds(tri.nfacets());

    for (auto h : tri.iter_halfedges()) {
		auto f = h.facet();
		auto opp_f = h.opposite().facet();

		if (tri_flag[f] != tri_flag[opp_f])
			continue;

        ds.merge(h.facet(), h.opposite().facet());
    }

    // Get associate facet id to group id
    std::vector<int> setIds;
    ds.get_sets_id(setIds);

	// Extract by groups
	std::map<int, std::vector<int>> element_by_group;
    for (long unsigned int i = 0; i < setIds.size(); i++) {
		element_by_group[setIds[i]].push_back(i);
    }

	std::vector<std::pair<int, UM::vec3>> flag_dirs;

	// Get bary of each group
	for (auto kv : element_by_group) {
		// Compute bary of all facets
		UM::vec3 bary{0,0,0};
		for (auto f_idx : kv.second) {
			Triangle3 t = UM::Surface::Facet(tri, f_idx);
			bary += t.bary_verts();
		}
		bary /= kv.second.size();
		std::cout << kv.first << ", bary: " << bary << std::endl;

		// Get flag of group
		int flag = tri_flag[kv.second.front()];
		flag_dirs.push_back({flag, bary});
	}

    std::cout << "n set: " << ds.nsets() << std::endl;
    std::cout << "n set: " << element_by_group.size() << std::endl;
	return flag_dirs;
}

void cell_facet_cross(UM::Hexahedra &hex, UM::Volume::Volume::Facet &start_f, std::function<void(UM::Volume::Facet&, int, int)> f) {
	assert(hex.connected());

	std::vector<bool> visited(24 * hex.ncells(), false);

	for (auto he : start_f.iter_halfedges()) {

		auto cur_he = he;

		int dist = 0;
		while (true) {

			auto opp_c = cur_he.opposite_f().opposite_c();
			if (!opp_c.active() || visited[cur_he])
				break;

			visited[cur_he] = true;

			auto n_he = opp_c.opposite_f().next().next();
			auto facet = n_he.facet();

			f(facet, he.id_in_facet(), dist++);
			cur_he = n_he;
		}
	}
}


// void loop_cut(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Facet&)> f) {
// 	assert(hex.connected());

// 	std::vector<bool> visited(24 * hex.ncells(), false);
// 	std::queue<int> q;
// 	q.push(start_he);
// 	visited[start_he] = true;


// 	while (!q.empty()) {

// 		auto he_idx = q.front();
// 		auto he = Volume::Halfedge(hex, he_idx);

// 		q.pop();

// 		auto opp_c = he.next().opposite_f().opposite_c();

// 		if (opp_c.active()) {
// 			auto n_he = opp_c.opposite_f().next();

// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		} else {
// 			auto n_he = he.next().opposite_f().next();
// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		}

// 		opp_c = he.prev().opposite_f().opposite_c();

// 		if (opp_c.active()) {
// 			auto n_he = opp_c.opposite_f().prev();

// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		} else {
// 			auto n_he = he.prev().opposite_f().prev();
// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		}

// 		// Propagate to front
// 		opp_c = he.opposite_f().next().next().opposite_f().opposite_c();

// 		if (opp_c.active()) {
// 			auto n_he = opp_c;

// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		} 

// 		// Process current
// 		auto facet = he.opposite_f().facet();
// 		f(facet);
		
// 	}
// }

// void loop_cut2(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge& /* cur_he */, bool/* on_border */)> f) {
// 	assert(hex.connected());
// 	std::vector<bool> visited(24 * hex.ncells(), false);

// 	auto cur_he = start_he;

// 	while (true) {

// 		f(cur_he, true);

// 		auto opp_c = cur_he.next().opposite_f().opposite_c();
// 		if (!opp_c.active() || cur_he.next().opposite_f().next().facet().on_boundary()) {
// 			cur_he = cur_he.next().opposite_f().next();

// 		} else {
// 			cur_he = opp_c.opposite_f().next();

// 			// on border ?
// 		}

// 		if (visited[cur_he])
// 			break;

// 		visited[cur_he] = true;


// 	}

// }

// void loop_cut2(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge& /* cur_he */, bool/* on_border */)> f) {
// 	assert(hex.connected());

// 	std::vector<bool> visited(24 * hex.ncells(), false);
// 	std::queue<int> q;
// 	q.push(start_he);
// 	visited[start_he] = true;

// 	while (!q.empty()) {

// 		auto he_idx = q.front();
// 		auto he = Volume::Halfedge(hex, he_idx);

// 		q.pop();

// 		bool on_boundary = !he.facet().on_boundary();

// 		for (auto around_he : he.opposite_f().facet().iter_halfedges()) {
// 			if (around_he.opposite_f().facet().on_boundary())
// 				f(around_he, on_boundary);
// 		}

// 		f(he, true);


// 		// dir: i
// 		auto opp_c = he.next().opposite_f().opposite_c();
// 		if (opp_c.active()) {
// 			auto n_he = opp_c.opposite_f().next();
// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		}
// 		// dir: -i
// 		opp_c = he.prev().opposite_f().opposite_c();
// 		if (opp_c.active()) {
// 			auto n_he = opp_c.opposite_f().prev();
// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		}

// 		// dir: j
// 		opp_c = he.opposite_f().next().next().opposite_f().opposite_c();
// 		if (opp_c.active()) {
// 			auto n_he = opp_c;
// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		}
// 		// dir: -j
// 		opp_c = he.opposite_c();
// 		if (opp_c.active()) {
// 			auto n_he = opp_c.opposite_f().next().next().opposite_f();
// 			if (!visited[n_he]) {
// 				q.push(n_he);
// 				visited[n_he] = true;
// 			}
// 		}

// 	}
// }


void loop_cut(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge&)> f) {
	assert(hex.connected());

	std::vector<bool> visited(24 * hex.ncells(), false);
	std::queue<int> q;
	q.push(start_he);
	visited[start_he] = true;

	while (!q.empty()) {

		auto he_idx = q.front();
		auto he = Volume::Halfedge(hex, he_idx);

		visited[he] = true;
		q.pop();

		auto opp_c = he.next().opposite_f().opposite_c();

		if (opp_c.active()) {
			auto n_he = opp_c.opposite_f().next();

			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		} else {
			auto n_he = he.next().opposite_f().next();
			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		}

		opp_c = he.prev().opposite_f().opposite_c();

		if (opp_c.active()) {
			auto n_he = opp_c.opposite_f().prev();

			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		} else {
			auto n_he = he.prev().opposite_f().prev();
			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		}

		// Process current
		f(he);
		
	}

}

void App::cursor_pos_callback(double x, double y, int source) {
	SimpleMeshApplicationExt::cursor_pos_callback(x, y, source);

	// Doesn't take account of hovering in camera mode
	if (context_.gui_mode == Camera || is_loading)
		return;

	// Try to pick cell
	context_.hovered_cell = pick(MESH_CELLS);

	// If a cell is hovered, try to pick cell edge / cell facet 
	// TODO mesh_metadata.cell_type == MESH_HEX is a quick fix, I should modfy pickup_cell_facet that just hold QUAD !
	if (context_.is_cell_hovered() && context_.mesh_metadata.cell_type == MESH_HEX) {
		index_t e_idx = pickup_cell_edge(picked_point_, context_.hovered_cell);
		auto [f_idx, lf_idx] = pickup_cell_facet(picked_point_, context_.hovered_cell);

		context_.hovered_edge = e_idx;
		context_.hovered_facet = f_idx;
		context_.hovered_lfacet = lf_idx;
	}

	if (context_.gui_mode == Painting && context_.left_mouse_pressed) {
		index_t f_idx = pick(MESH_FACETS);
		// auto [cf_idx, _] = pickup_cell_facet2(context_.click_pos, context_.hovered_cell);

		if (f_idx != NO_FACET && f_idx < mesh_.facets.nb()) {

			GEO::Attribute<GEO::signed_index_t> flag(
				mesh_.facets.attributes(), "flag"
			);

			flag[f_idx] = context_.paint_value;

		}
	}
	else if (context_.gui_mode == LayerPadding) {
		layer_pad_tool.hover_callback(x, y, source);
	}
	else if (context_.gui_mode == BlocPadding) {
		bloc_pad_tool.hover_callback(x, y, source);
	}
}

void App::mouse_button_callback(int button, int action, int mods, int source) {

	if (context_.gui_mode == Camera) {
    	SimpleMeshApplication::mouse_button_callback(button,action,mods,source);
		return;
	}

	if (action == EVENT_ACTION_DOWN && button == 0) {
		context_.left_mouse_pressed = true;
		context_.click_pos = picked_point_;
	}
	else if (action == EVENT_ACTION_UP && button == 0) {
		context_.left_mouse_pressed = false;
	}


	if (action == EVENT_ACTION_DOWN && button == 0) {
		context_.selected_vertex = context_.hovered_vertex;
		context_.selected_edge = context_.hovered_edge;
		context_.selected_facet = context_.hovered_facet;
		context_.selected_lfacet = context_.hovered_lfacet;
		context_.selected_cell = context_.hovered_cell;
	}

	// If left click
    if (context_.gui_mode == Painting && action == EVENT_ACTION_DOWN && button == 0) {


    }
	else if (context_.gui_mode == BlocPadding && action == EVENT_ACTION_UP && button == 0) {
		bloc_pad_tool.mouse_button_callback(button, action, mods, source);
	}
	else if (context_.gui_mode == Painting && action == EVENT_ACTION_UP && button == 0) {
		// Transfert attribute from surface tri to volume tet
		FacetAttribute<int> tri_flag(context_.tet_bound.tri, -1);
		CellFacetAttribute<int> tet_flag(context_.tet_bound.tet, -1);
		um_bindings::um_attr_from_geo_attr<GEO::MESH_FACETS>(mesh_, "flag", context_.tet_bound.tri, tri_flag.ptr);

		context_.tet_bound.set_attribute_to_volume(tri_flag, tet_flag);
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(context_.tet_bound.tet, tet_flag.ptr, "tet_flag", mesh_);
	}
	else if (context_.gui_mode == LayerPadding && action == EVENT_ACTION_UP && button == 0) {
		layer_pad_tool.mouse_button_callback(button, action, mods, source);
	}



}

void App::key_callback(int key, int scancode, int action, int mods) {
	
	std::cout << "key pressed: " << key << std::endl;

	// Ctrl
	if (action == EVENT_ACTION_DOWN && key == 341) {
		context_.switch_mode = context_.gui_mode;
		context_.gui_mode = Camera;
	} else if (action == EVENT_ACTION_UP && key == 341) {
		context_.gui_mode = context_.switch_mode;
		context_.switch_mode = Camera;
	}

	// Validate loop pad
	if (context_.gui_mode == LayerPadding && (key == 257 || key == 335) && action == EVENT_ACTION_DOWN) {
		layer_pad_tool.validate_callback();
	}
	else if (context_.gui_mode == BlocPadding && (key == 257 || key == 335) && action == EVENT_ACTION_DOWN) {
		bloc_pad_tool.validate_callback();
	}
	else if (context_.gui_mode == BlocPadding && (key == 256) && action == EVENT_ACTION_DOWN) {
		bloc_pad_tool.escape_callback();
	}


}

bool App::save(const std::string& filename) {
    // if(String::string_ends_with(filename,".txt")) { // bypass inherited save behavior in case of a .txt file -> save the labeling only
    //     save_labeling(filename,mesh_,labeling_);
    //     fmt::println(Logger::out("I/O"),"Labeling saved to {}",filename); Logger::out("I/O").flush();
    //     return true;
    // }
	return SimpleMeshApplication::save(filename);
}

void App::reset() {
	// Reset all values
	context_.gui_mode = Camera;
	// Clear selections
	context_.reset_hovered_selected();

	// Clear tools
	layer_pad_tool.clear();
	bloc_pad_tool.clear();

	// Clear UM meshes
	context_.tet.clear();
	context_.hex.clear();
}

bool App::load(const std::string& filename) {
    
	mesh_gfx_.set_mesh(nullptr);
    mesh_.clear(false, true);

    MeshIOFlags flags;
	
	is_loading = true;

	std::string mesh_filename = filename;
	std::filesystem::path filename_path(filename);


	// Load mesh metadata from json file ?
	if (filename_path.extension() == ".json") {
		std::ifstream ifs(filename);
		std::string content;
		ifs >> content;
		ifs.close();

		auto json = json::parse(content);
		context_.mesh_metadata = MeshMetadata::from_json(json);
		mesh_filename = context_.mesh_metadata.filename;
	} else {
		context_.mesh_metadata.filename = filename;
		context_.mesh_metadata.cell_type = MESH_TET;
	}

    if(!mesh_load(mesh_filename, mesh_, flags)) {
		// TODO clear all here too !
        return false;
    }

	reset();



	// Init UM tet from GEO mesh
	if (context_.mesh_metadata.cell_type == GEO::MESH_TET) {
		um_bindings::um_tet_from_geo_mesh(mesh_, context_.tet);
		context_.tet.connect();
	}
	else if (context_.mesh_metadata.cell_type == GEO::MESH_HEX) {
		um_bindings::um_hex_from_geo_mesh(mesh_, context_.hex);
		context_.hex.connect();
	}

	is_loading = false;


	// Display info
    mesh_gfx_.set_animate(false);
    mesh_.vertices.set_dimension(3);
    mesh_.vertices.set_double_precision(); // added
    mesh_gfx_.set_mesh(&mesh_);
    current_file_ = mesh_filename;

    labeling_visu_mode_transition();

    clear_scene_overlay();

	std::cout << "Mesh loaded !" << std::endl;

    return true;
}


void App::draw_object_properties() {

	ImGui::Checkbox("Tool preview", &tool_preview);

	int n_facet_per_cell = context_.mesh_metadata.cell_type == MESH_HEX ? 6 : 4;

	ImGui::Text("Hovered vertex: %i", context_.hovered_vertex);
	ImGui::Text("Hovered cell edge: %i", context_.hovered_edge);
	ImGui::Text("Hovered cell facet: %i", context_.hovered_facet);
	ImGui::Text("Hovered cell facet [UM]: %i", um_bindings::um_facet_index_from_geo_facet_index(context_.hovered_facet, n_facet_per_cell));
	ImGui::Text("Hovered cell local facet: %i", context_.hovered_lfacet);
	ImGui::Text("Hovered cell: %i", context_.hovered_cell);

	// ImGui::Text("Hen: %i", he_n);
	// if(ImGui::Button("Next!")) {
	// 	he_n++;
	// }

	if(ImGui::Button("Camera")) {
		context_.gui_mode = Camera;
	}

	if(ImGui::Button("Hover")) {
		context_.gui_mode = Hover;
	}

	ImGui::Separator();

	bool is_visible_compute_flag_tool = (context_.mesh_metadata.cell_type == MESH_TET);

	if (is_visible_compute_flag_tool) {


		// auto t = convert_to_ImTextureID(colormaps_[current_colormap_index_].texture);
		// ImGui::ColorButton("+X", ImVec4(1,0,0,1));
		ImGui::TextUnformatted("Paint flags");
		if(ImGui::Button("No Paint")) {
			context_.paint_value = -1;
			context_.gui_mode = Painting;
		}
		if(ImGui::Button("-X")) {
			context_.paint_value = 0;
			context_.gui_mode = Painting;
		}
		ImGui::SameLine();
		if(ImGui::Button("-Y")) {
			context_.paint_value = 1;
			context_.gui_mode = Painting;
		}
		ImGui::SameLine();
		if(ImGui::Button("-Z")) {
			context_.paint_value = 2;
			context_.gui_mode = Painting;
		}	
		if(ImGui::Button("+X")) {
			context_.paint_value = 3;
			context_.gui_mode = Painting;
		}
		ImGui::SameLine();
		if(ImGui::Button("+Y")) {
			context_.paint_value = 4;
			context_.gui_mode = Painting;
		}
		ImGui::SameLine();
		if(ImGui::Button("+Z")) {
			context_.paint_value = 5;
			context_.gui_mode = Painting;
		}	

		ImGui::Separator();

	}


	if (layer_pad_tool.is_compatible()) {
		if (layer_pad_tool.draw_gui()) {
			// Clean other tools
			bloc_pad_tool.clear();
		}
	}

	if (bloc_pad_tool.is_compatible()) {
		if (bloc_pad_tool.draw_gui()) {
			// Clean other tools
			layer_pad_tool.clear();
		}
	}	

	if (is_visible_compute_flag_tool) {

		if (ImGui::Button("Compute flags !")) {
			
			// Compute flag on tet and tri
			// TetBoundary tet_bound(tet);
			context_.tet_bound.update();
			
			// To GEO mesh
			// TODO maybe move under tet_bound.set_attribute_to_surface(tet_flag, tri_flag);
			um_bindings::geo_mesh_from_tetboundary(context_.tet_bound, mesh_);

			UM::CellFacetAttribute<int> tet_flag(context_.tet, -1);
			UM::FacetAttribute<int> tri_flag(context_.tet_bound.tri, -1);

			// Compute flag
			algo::naive_tag(context_.tet, tet_flag);
			// Transfert flag from tet to tri for display
			context_.tet_bound.set_attribute_to_surface(tet_flag, tri_flag);
			// Update GEO mesh attribute "flag"
			um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(context_.tet, tet_flag.ptr, "tet_flag", mesh_);
			um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(context_.tet_bound.tri, tri_flag.ptr, "flag", mesh_);


			// TODO encapsulate in atomic unit ! + try catch to guarentee consistency
			// Save mesh metadata
			context_.mesh_metadata = { 
				.filename = "flagged.geogram", 
				.cell_type = GEO::MESH_TET, 
				.attributes = {
					{
						.name = "tet_flag", 
						.type = "int", 
						.where = MESH_CELL_FACETS
					}
				} 
			};
			// Write mesh
			write_by_extension(context_.mesh_metadata.filename, context_.tet_bound.tet, {{}, {}, {{"tet_flag", tet_flag.ptr}}, {}});
			context_.mesh_metadata.save();

			labeling_visu_mode_transition();
			show_surface_ = true;
			show_volume_ = false;
		}

		if (ImGui::Button("Compute patches !")) {
			// Compute flag on tet and tri
			TetBoundary tet_bound(context_.tet);
			UM::FacetAttribute<int> tri_flag(tet_bound.tri, -1);
			um_bindings::um_attr_from_geo_attr<GEO::MESH_FACETS>(mesh_, "flag", tet_bound.tri, tri_flag.ptr);
			flag_dirs = compute_patches(tet_bound.tri, tri_flag);
		}

		ImGui::Separator();

	}


	// Criteria to display polycubify tool
	auto mesh_metadata_attr = context_.mesh_metadata.get_attr("tet_flag");
	
	bool is_visible_polycubify_tool = (context_.mesh_metadata.cell_type == MESH_TET) && (mesh_metadata_attr.has_value() && mesh_metadata_attr.value().where == GEO::MESH_CELL_FACETS);
	
	if (is_visible_polycubify_tool) {
		int nhex_wanted = 3000;
		if (ImGui::Button("Polycubify !")) {

			// Get UM cell facet attribute tet_flag from GEO mesh
			UM::CellFacetAttribute<int> tet_flag(context_.tet, -1);
			um_bindings::um_attr_from_geo_attr<GEO::MESH_CELL_FACETS>(mesh_, "tet_flag", context_.tet, tet_flag.ptr);

			try {
				BenjaminAPI::polycubify(context_.tet, tet_flag, context_.hex, nhex_wanted);
			} catch (const std::runtime_error &e) {
				Logger::warn("An error occur when trying to polycubify. Detail: " + std::string(e.what()));
				std::cout << "polycubify fail" << std::endl;
				return;
			}

			HexBoundary hex_bound(context_.hex);
			// Replace current GEO mesh by UM Hex
			um_bindings::geo_mesh_from_hexboundary(hex_bound, mesh_);


			// write_by_extension("input.geogram", tet);
			// write_by_extension("polycubify_hex.geogram", hex);
			// write_by_extension("polycubify_quad.geogram", hex_bound.quad);
			// mesh_save(mesh_, "polycubify2.geogram");
			// TODO encapsulate in atomic unit ! + try catch to guarentee consistency

			// Write mesh

			// Save mesh metadata in json !!!!
			context_.mesh_metadata = { 
				.filename = "polycubified.geogram", 
				.cell_type = GEO::MESH_HEX, 
				.attributes = {} 
			};
			write_by_extension(context_.mesh_metadata.filename, hex_bound.hex, {{}, {}, {}, {}});
			context_.mesh_metadata.save();

			// View
			mesh_gfx_.set_mesh(&mesh_);
			labeling_visu_mode_transition();
			show_surface_ = false;
			show_volume_ = true;

		}

	}

}

std::string App::supported_write_file_extensions() {
    return SimpleMeshApplication::supported_write_file_extensions() + ";json"; // add .json in supported write file extensions
}

std::string App::supported_read_file_extensions() {
    return SimpleMeshApplication::supported_read_file_extensions() + ";json"; // add .json in supported write file extensions
}


void App::labeling_visu_mode_transition() {
    if(colormaps_.empty()) {
        // GL is not initialized yet
        // the state transition will be triggered later, in GL_initialize()
        std::cout << "empty colormap ! " << std::endl;
        return;
    }

	// TODO Add view mode for mesh type to switch view

	show_attributes_ = true;
    show_vertices_ = false;
	show_volume_ = true;
	show_surface_ = true;
	show_hexes_ = true;
	
	current_colormap_index_ = COLORMAP_RAINBOW;
	attribute_ = "facets.flag";
	attribute_subelements_ = MESH_FACETS;
	attribute_name_ = "flag";
	attribute_min_ = -1;
	attribute_max_ = 5;
	lighting_ = false;
}