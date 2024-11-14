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


using json = nlohmann::json;

App::App(const std::string name) : 
	SimpleMeshApplicationExt(name), 
	context_(mesh_, mesh_gfx_, {
		.show_attributes_ = show_attributes_,
		.show_vertices_ = show_vertices_,
		.show_volume_ = show_volume_,
		.show_surface_ = show_surface_,
		.show_hexes_ = show_hexes_,
		.current_colormap_index_ = current_colormap_index_,
		.attribute_ = attribute_,
		.attribute_subelements_ = attribute_subelements_,
		.attribute_name_ = attribute_name_,
		.attribute_min_ = attribute_min_,
		.attribute_max_ = attribute_max_,
		.lighting_ = lighting_
	}),
	hover_selection_colors_({
		{1.f,0.3f,0.6f, 1.f}, // pink hover
		{0.95f,0.2f,0.0f, 1.f}  // red select
	}),
	// 255 128 128
	//
	flagging_colors_({
		{0.5f, 0.5f, 0.5f, 1.0f},
		{1.f, 0.6f, 0.6f, 1.0f},
		{0.6f, 1.0f, 0.6f, 1.0f},
		{0.6f, 0.6f, 1.0f, 1.0f},
		{0.8f, 0.16f, 0.16f, 1.0f},
		{0.16f, 0.8f, 0.16f, 1.0f},
		{0.16f, 0.16f, 0.8f, 1.0f}
	}),
	camera_tool(context_),
	hover_tool(context_),
	layer_pad_tool(context_),
	bloc_pad_tool(context_),
	new_bloc_pad_tool(context_),
	paint_flag_tool(context_),
	polycubify_tool(context_)
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

void App::normalize_mesh() {
	double xyz_min[3];
	double xyz_max[3];
	get_bbox(mesh_, xyz_min, xyz_max, false);
	GEO::vec3 bb_min(xyz_min[0], xyz_min[1], xyz_min[2]);
	GEO::vec3 bb_max(xyz_max[0], xyz_max[1], xyz_max[2]);
	double max_coord = std::max(std::max(xyz_max[1], xyz_max[2]), xyz_max[0]);
	std::cout << "min: " << bb_min << ", max: " << bb_max << std::endl;
	std::cout << "max coord: " << max_coord << std::endl;

	// normalize
	for (int v = 0; v < mesh_.vertices.nb(); v++) {
		mesh_.vertices.point(v) /= max_coord; 
		// mesh_.vertices.point(v) *= 2; 
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
	tools[context_.gui_mode]->draw(hovered_color, selected_color, hover_selection_colormap);
	// if (context_.gui_mode == Hover) {
	// 	hover_tool.draw(hovered_color, selected_color, hover_selection_colormap);
	// }
	// if (context_.gui_mode == LayerPadding) {
	// 	layer_pad_tool.draw(hovered_color, selected_color, hover_selection_colormap);
	// } 
	// if (context_.gui_mode == BlocPadding) {
	// 	bloc_pad_tool.draw(hovered_color, selected_color, hover_selection_colormap);
	// }



	// Last picked point position as point
	if (context_.show_last_picked_point_) {
		glupBegin(GLUP_POINTS);
		glupPrivateVertex3dv(context_.click_pos.data());
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

	if (show_features) {
		glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(0., 0., 1., 1.).data());
		glupBegin(GLUP_LINES);
		for (int e = 0; e < mesh_.edges.nb(); e++) {
			auto p0 = mesh_.vertices.point(mesh_.edges.vertex(e, 0));
			auto p1 = mesh_.vertices.point(mesh_.edges.vertex(e, 1));
			

			glupPrivateVertex3dv(p0.data());
			glupPrivateVertex3dv(p1.data());
		}
		glupEnd();
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

	ImGui::InputText("Gmhsh path", gmsh_path, sizeof(gmsh_path));

	if (load_step) {
		ImGui::OpenPopup("Load step options");
		if (ImGui::BeginPopupModal("Load step options", NULL, ImGuiWindowFlags_AlwaysAutoResize))
		{

			// bool volume;
			
			ImGui::InputFloat("Mesh size factor", &size_factor);
			// ImGui::Checkbox("tetrahedrize ?", &volume);

			if (ImGui::Button("Load", ImVec2(120, 0)))
			{
				std::string out_filename = context_.mesh_metadata.filename + ".mesh";
				// Write config file for GMSH
				std::ofstream conf_file;
				conf_file.open("conf.geo");
				if (conf_file.is_open()) { 
					// std::string mesh_type = volume ? "3" : "2";
					std::string mesh_type = "3";
					conf_file << "General.Terminal = 1;\nMesh.MeshSizeFactor = " + std::to_string(size_factor) + ";\nMesh.AngleToleranceFacetOverlap = 0.02;\nMesh " + mesh_type + ";\nSave \"" + out_filename + "\";";
					conf_file.close();
				} else {
					std::cout << "Unable to write file: conf.geo at the Graphite root directory." << std::endl;
					return;
				}

				std::string command = std::string(gmsh_path) + " " + context_.mesh_metadata.filename + " conf.geo -0";
				int result = system(command.c_str());

				if (result != 0) {
					std::cout << "Unable to mesh step file from GMSH." << std::endl;
					return;
				}

				// TODO refactor this
				MeshIOFlags flags;
				if(!mesh_load(out_filename, mesh_, flags)) {
					reset();
					return;
				}

				reset();

				normalize_mesh();

				// Init UM tet from GEO mesh
				context_.mesh_metadata.cell_type = MESH_TET;
				um_bindings::um_tet_from_geo_mesh(mesh_, context_.tet);
				context_.tet.connect();

				is_loading = false;


				// Display info
				mesh_gfx_.set_animate(false);
				mesh_.vertices.set_dimension(3);
				mesh_.vertices.set_double_precision(); // added
				mesh_gfx_.set_mesh(&mesh_);
				current_file_ = out_filename;

				labeling_visu_mode_transition();

				clear_scene_overlay();

				std::cout << "Mesh loaded !" << std::endl;

				load_step = false;
				ImGui::CloseCurrentPopup();
			}

			ImGui::EndPopup();
		}
	}



	ImGui::Checkbox("Show grid", &show_grid_);
	ImGui::Checkbox("Show axes", &show_axes_);
	ImGui::Checkbox("Show picked point", &context_.show_last_picked_point_);
	ImGui::Separator();
	ImGui::Checkbox("Show vertices", &show_vertices_);
	ImGui::Checkbox("Show features", &show_features);

	// ImGui::Checkbox("Show surface", &show_surface_);
	// ImGui::Checkbox("Show volume", &show_volume_);
	ImGui::TextUnformatted("View mode");
	for (int i = 0; i < IM_ARRAYSIZE(context_.view.modes); i++)
	{
		bool isSelected = (context_.view.current_mode == i);
		if (ImGui::Selectable(context_.view.modes[i], isSelected))
		{
			context_.view.change_mode(i);
		}

		if (isSelected)
			ImGui::SetItemDefaultFocus();
	}

	ImGui::Separator();

	tools[context_.gui_mode]->draw_viewer_properties();
}

void App::GL_initialize() {
    SimpleMeshApplicationExt::GL_initialize();
	init_rgba_colormap("hover_selection", 2, 1, hover_selection_colors_.as_chars());
	init_rgba_colormap("flagging", 7, 1, flagging_colors_.as_chars());
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



void App::cursor_pos_callback(double x, double y, int source) {
	SimpleMeshApplicationExt::cursor_pos_callback(x, y, source);

	// Doesn't take account of hovering in camera mode
	if (context_.gui_mode == Camera || is_loading)
		return;

	if ((UM::vec2{x,y} - last_mouse_pos).norm2() < 9)
		return;

	// Try to pick cell
	context_.hovered_cell = pick(MESH_CELLS);

	// Get facet idx only when surface is displayed
	if (context_.view.current_mode == ViewBinding::Mode::Surface)
		context_.hovered_facet = pick(MESH_FACETS);

	// If a cell is hovered, try to pick cell edge / cell facet 
	if (context_.view.current_mode == ViewBinding::Mode::Volume && context_.is_cell_hovered()) {
		index_t e_idx = pickup_cell_edge(picked_point_, context_.hovered_cell);
		context_.hovered_edge = e_idx;

		auto [f_idx, lf_idx] = pickup_cell_facet(picked_point_, context_.hovered_cell);
		context_.hovered_cell_facet = f_idx;
		context_.hovered_cell_lfacet = lf_idx;
	}

	tools[context_.gui_mode]->hover_callback(x, y, source);

	last_mouse_pos = {x, y};
}

void App::mouse_button_callback(int button, int action, int mods, int source) {
	// std::cout << "mouse: " << button << "," << action << std::endl;

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
		// TODO add function select_hovered
		context_.selected_vertex = context_.hovered_vertex;
		context_.selected_edge = context_.hovered_edge;
		context_.selected_facet = context_.hovered_facet;
		context_.selected_cell_facet = context_.hovered_cell_facet;
		context_.selected_cell_lfacet = context_.hovered_cell_lfacet;
		context_.selected_cell = context_.hovered_cell;
	}

	if (action == EVENT_ACTION_UP && button == 0) {
		tools[context_.gui_mode]->mouse_button_callback(button, action, mods, source);
	}

}

void App::scroll_callback(double xoffset, double yoffset) {
	if (context_.gui_mode == Camera) {
    	SimpleMeshApplication::scroll_callback(xoffset, yoffset);
		return;
	}

	tools[context_.gui_mode]->scroll_callback(xoffset, yoffset);
}

void App::key_callback(int key, int scancode, int action, int mods) {
	
	// std::cout << "key pressed: " << key << std::endl;

	// Ctrl
	if (action == EVENT_ACTION_DOWN && key == 341) {
		context_.switch_mode = context_.gui_mode;
		context_.gui_mode = Camera;
	} else if (action == EVENT_ACTION_UP && key == 341) {
		context_.gui_mode = context_.switch_mode;
		context_.switch_mode = Camera;
	}

	// Validate
	if ((key == 257 || key == 335) && action == EVENT_ACTION_DOWN) {
		std::cout << "gui mode: " << context_.gui_mode << std::endl;
		tools[context_.gui_mode]->validate_callback();
	} else if (key == 256 && action == EVENT_ACTION_DOWN) {
		tools[context_.gui_mode]->escape_callback();
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
	for (auto &tool : tools) {
		tool->clear();
	}

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


	// Load mesh metadata from json file !
	if (filename_path.extension() == ".json") {
		std::ifstream ifs(filename);
		std::string content;
		ifs >> content;
		ifs.close();

		auto json = json::parse(content);
		context_.mesh_metadata = MeshMetadata::from_json(json);
		mesh_filename = context_.mesh_metadata.filename;
	} else if (filename_path.extension() == ".step" || filename_path.extension() == ".stp") {

		context_.mesh_metadata.filename = filename;
		load_step = true;
		return true;


	} else {
		context_.mesh_metadata.filename = filename;
		context_.mesh_metadata.cell_type = MESH_TET;
	}

    if(!mesh_load(mesh_filename, mesh_, flags)) {
		reset();
        return false;
    }
	
	reset();

	normalize_mesh();

	// Init UM tet from GEO mesh
	if (context_.mesh_metadata.cell_type == GEO::MESH_TET) {
		um_bindings::um_tet_from_geo_mesh(mesh_, context_.tet);
		context_.tet.connect();
		context_.view.change_mode(ViewBinding::Mode::Surface);
	}
	else if (context_.mesh_metadata.cell_type == GEO::MESH_HEX) {
		um_bindings::um_hex_from_geo_mesh(mesh_, context_.hex);
		context_.hex.connect();
		context_.view.change_mode(ViewBinding::Mode::Volume);
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

	ImGui::Text("Selected tool: %s", tools[context_.gui_mode]->get_name().c_str());
	ImGui::Text("Hovered vertex: %i", context_.hovered_vertex);
	ImGui::Text("Hovered cell edge: %i", context_.hovered_edge);
	ImGui::Text("Hovered facet: %i", context_.hovered_facet);
	ImGui::Text("Hovered cell facet: %i", context_.hovered_cell_facet);
	ImGui::Text("Hovered cell facet [UM]: %i", um_bindings::um_facet_index_from_geo_facet_index(context_.hovered_cell_facet, n_facet_per_cell));
	ImGui::Text("Hovered cell local facet: %i", context_.hovered_cell_lfacet);
	ImGui::Text("Hovered cell: %i", context_.hovered_cell);

	// ImGui::Text("Hen: %i", he_n);
	// if(ImGui::Button("Next!")) {
	// 	he_n++;
	// }

	for (auto &tool : tools) {
		if (tool->is_compatible()) {
			tool->draw_object_properties();
		}
	}

}

std::string App::supported_write_file_extensions() {
    return SimpleMeshApplication::supported_write_file_extensions() + ";json"; // add .json in supported write file extensions
}

std::string App::supported_read_file_extensions() {
    return SimpleMeshApplication::supported_read_file_extensions() + ";json;step;stp"; // add .json, .step in supported write file extensions
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
	// show_volume_ = true;
	// show_surface_ = true;
	show_hexes_ = true;
	
	current_colormap_index_ = COLORMAP_FLAGGING;
	attribute_ = "facets.flag";
	attribute_subelements_ = MESH_FACETS;
	attribute_name_ = "flag";
	attribute_min_ = -1;
	attribute_max_ = 5;
	lighting_ = false;
}