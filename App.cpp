#include <ultimaille/all.h>
#include <geogram_gfx/full_screen_effects/ambient_occlusion.h>

#include "App.h"

#include "geom_ultimaille_binding.h"
#include "tag_face.h"

#include <nicostuff/algo/framework/benjamin_API.h>
#include "my_tet_boundary.h"

#include <queue>

App::App(const std::string name) : SimpleMeshApplicationExt(name)/*, mesh_binding(mesh_)*/ {

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

inline Quaternion align_with_uv(UM::vec3 u, UM::vec3 v) {
	v.normalize();
	u.normalize();
	
	if (u * v < -.99999) { 
		Quaternion res;  
		res.v = {1.0, 0., 0.}; 
		res.w = 0;
		return res; 
	}
	
	if (std::abs(u * v) > .99999)
		return Quaternion();

	UM::vec3 inbetwen_uv(v + u);
	inbetwen_uv.normalize();

	Quaternion res;
	res.w = v * inbetwen_uv; // scalar product with (1,0,0) divided by norm
	res.v = cross(inbetwen_uv, v); // cross product with (1,0,0) 
	
	return res;
}

void draw_arrow(UM::vec3 A, UM::vec3 B, double r, int res, double cone_tube_ratio, GEO::vec4f color) {
	
	if (r == -1) 
		r = (B - A).norm() / 10.;

	UM::vec3 n = B - A;
	double l = std::max(1e-5, n.norm());
	n = n / l;

	Quaternion quat = align_with_uv({0, 0, 1}, n);
	auto M = quat.rotation_matrix();

	glupSetColor4fv(GLUP_FRONT_COLOR, color.data());
	auto mesh_width = glupGetMeshWidth();
	glupSetMeshWidth(0.);

	// Set temporary picking mode to constant to deactivate picking on GL primitive
	auto picking_mode = glupGetPickingMode();
	glupPickingMode(GLUP_PICK_CONSTANT);

	glupBegin(GLUP_TRIANGLES);

	for (int v = 0; v < res; v++) {
		// circle parametric equation
		auto pt = r * UM::vec3{cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0.};
		auto pt_nxt = r * UM::vec3{cos(2. * M_PI * double(v + 1) / double(res)), sin(2. * M_PI * double(v + 1) / double(res)), 0.};

		auto pt_base = (r * cone_tube_ratio) * UM::vec3{cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0.};
		auto pt_nxt_base = (r * cone_tube_ratio) * UM::vec3{cos(2. * M_PI * double(v + 1) / double(res)), sin(2. * M_PI * double(v + 1) / double(res)), 0.};

		auto apex = UM::vec3{0, 0, l};

		auto h = UM::vec3{0, 0, std::max(0., l - 2.5 * r)};
		auto pt_cone_base = pt + h;
		auto pt_nxt_cone_base = pt_nxt + h;

		// Cone
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_cone_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_cone_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * apex + A).data());

		// Cone base
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_cone_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_cone_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * h + A).data());

		// Tube base
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * UM::vec3{0,0,0} + A).data());

		// Tube
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_base + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * (pt_nxt_base + h) + A).data());

		glupPrivateVertex3dv(um_bindings::geo_vec(M * (pt_nxt_base + h) + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * (pt_base + h) + A).data());
		glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_base + A).data());

	}

	glupEnd();
	glupSetMeshWidth(mesh_width);
	glupPickingMode(picking_mode);
}

void App::draw_path(std::vector<UM::vec3> &path, GEO::vec4f color, bool arrow) {

	auto picking_mode = glupGetPickingMode();
	glupPickingMode(GLUP_PICK_CONSTANT);

	glupSetColor4fv(GLUP_FRONT_COLOR, color.data());
	glupBegin(GLUP_POINTS);
	for (auto p : path) {
		double pp[3] = {p.x, p.y, p.z};
		glupPrivateVertex3dv(pp);
	}
	glupEnd();

	glupPickingMode(picking_mode);

	if (!arrow) {
		glupBegin(GLUP_LINES);
		for (int i = 0; i < path.size(); i+=2) {
			auto pA = path[i];
			auto pB = path[i + 1];
			double ppA[3] = {pA.x, pA.y, pA.z};
			double ppB[3] = {pB.x, pB.y, pB.z};
			glupPrivateVertex3dv(ppA);
			glupPrivateVertex3dv(ppB);
		}
		glupEnd();
	} else {
		for (int i = 0; i < path.size(); i+=2) {
			draw_arrow(path[i], path[i + 1], -1, 10, 0.5, color);
		}
	}

}

void App::draw_axis() {
	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(1.f,0.f,0.f).data());
	glupBegin(GLUP_LINES);
	glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
	glupPrivateVertex3dv(GEO::vec3(1, 0, 0).data());
	glupEnd();

	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(0.f,1.f,0.f).data());
	glupBegin(GLUP_LINES);
	glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
	glupPrivateVertex3dv(GEO::vec3(0, 1, 0).data());
	glupEnd();

	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(0.f,0.f,1.f).data());
	glupBegin(GLUP_LINES);
	glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
	glupPrivateVertex3dv(GEO::vec3(0, 0, 1).data());
	glupEnd();

	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(1.f, 1.f, 1.f).data());
	glupBegin(GLUP_POINTS);
	glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
	glupEnd();
}

void App::draw_grid() {
	double scale = 10.;
	int n_lines = 9;
	double y = -.1;

	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(0.2f,0.2f,0.2f).data());
	glupBegin(GLUP_LINES);
	for (int i = 0; i < n_lines; i++) {
		double j = i * (scale / n_lines) - scale * .5;
		double k = scale * .5;

		GEO::vec3 a(j, y, k);
		GEO::vec3 b(j, y, -k);
		
		glupPrivateVertex3dv(a.data());
		glupPrivateVertex3dv(b.data());
	}
	for (int i = 0; i < n_lines; i++) {
		double j = i * (scale / n_lines) - scale * .5;
		double k = scale * .5;

		GEO::vec3 a(k, y, j);
		GEO::vec3 b(-k, y, j);
		
		glupPrivateVertex3dv(a.data());
		glupPrivateVertex3dv(b.data());
	}
	glupEnd();
}

void App::draw_scene() {
    SimpleMeshApplicationExt::draw_scene();

	// GLUPfloat *view;
	// glupGetMatrixfv(GLUP_MODELVIEW_MATRIX, view);
	// ImGuizmo::ViewManipulate(view, 2.f, ImVec2(100,100), ImVec2(20,20), ImU32(100));

	glupSetPointSize(10.0);

	if (display_axes_)
		draw_axis();
	if (display_grid_)
		draw_grid();

	

	// Display hovered halfedge

	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(1.f,0.3f,0.6f, 1.f).data());
	glupBegin(GLUP_LINES);
	glupPrivateVertex3dv(posA.data());
	glupPrivateVertex3dv(posB.data());
	glupEnd();

	glupBegin(GLUP_POINTS);
	glupPrivateVertex3dv(posA.data());
	glupPrivateVertex3dv(posB.data());
	glupEnd();

	// PATH
	draw_path(hovered_path, GEO::vec4f(1.f,0.3f,0.6f, 1.f), true);
	draw_path(selected_path, GEO::vec4f(1.f,0.2f,0.0f, 1.f), false);


	// Display UM halfedges to remove !
	// glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(1.f,0.0f,0.0f, 1.f).data());
	// glupBegin(GLUP_POINTS);
	// double ppA[3] = {posAb.x, posAb.y, posAb.z};
	// glupPrivateVertex3dv(ppA);
	// glupEnd();

	// glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(0.f,1.0f,0.0f, 1.f).data());
	// glupBegin(GLUP_POINTS);
	// double ppB[3] = {posBb.x, posBb.y, posBb.z};
	// glupPrivateVertex3dv(ppB);
	// glupEnd();

	// glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(0.f,1.0f,0.0f, 1.f).data());
	// glupBegin(GLUP_LINES);
	// auto mid = (posAb + posBb) * .5;
	// auto n = mid + posN;
	// double pMid[3] = {mid.x, mid.y, mid.z};
	// double pN[3] = {n.x, n.y, n.z};
	// glupPrivateVertex3dv(pMid);
	// glupPrivateVertex3dv(pN);
	// glupEnd();

	glupBegin(GLUP_POINTS);
	glupPrivateVertex3dv(click_pos.data());
	glupEnd();



	// if (hovered_facet != NO_FACET && hovered_facet < mesh_.cell_facets.nb()) {
	// 	glupSetPointSize(5.0);
	// 	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(0.f, 0.2f, 0.6f, 0.2f).data());
	// 	glupBegin(GLUP_POINTS);
		
	// 	auto facet_n_verts = mesh_.cells.facet_nb_vertices(hovered_cell, hovered_lfacet);

	// 	for (int lv = 0; lv < facet_n_verts; lv++) {
	// 		auto a = mesh_.cells.facet_vertex(hovered_cell, hovered_lfacet, lv);
	// 		glupPrivateVertex3dv(mesh_.vertices.point(a).data());
	// 	}
	// 	// for (int lv = 0; lv < facet_n_verts; lv++) {
			
	// 	// 	auto a = mesh_.cells.facet_vertex(hovered_cell, hovered_lfacet, lv);
	// 	// 	auto b = mesh_.cells.facet_vertex(hovered_cell, hovered_lfacet, (lv + 1) % facet_n_verts);
	// 	// 	auto c = mesh_.cells.facet_vertex(hovered_cell, hovered_lfacet, (lv + 2) % facet_n_verts);
	// 	// 	glupPrivateVertex3dv(mesh_.vertices.point(a).data());
	// 	// 	glupPrivateVertex3dv(mesh_.vertices.point(b).data());
	// 	// 	glupPrivateVertex3dv(mesh_.vertices.point(c).data());
	// 	// }
	// 	glupEnd();
	// }


	// EXAMPLE TRIANGLES
	// glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(rand()% 100 / 100.,rand()% 100 / 100.,rand()% 100 / 100.).data());
	// glupBegin(GLUP_TRIANGLES);
	// for (int i = 0; i < 20; i++) {
	// 	glupPrivateVertex3dv(GEO::vec3(i / 10., 0, 0).data());
	// 	glupPrivateVertex3dv(GEO::vec3(i / 10., 1, 0).data());
	// 	glupPrivateVertex3dv(GEO::vec3(i / 10., 1, 1).data());
	// }
	// glupEnd();



	// EXAMPLE LINES TEXTURED
	// if (tet.connected()) {
	// 	glupEnable(GLUP_TEXTURING);
	// 	glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
	// 	glBindTexture(GL_TEXTURE_2D, colormaps_[COLORMAP_PARULA].texture);
	// 	glupTextureType(GLUP_TEXTURE_2D);
	// 	glupTextureMode(GLUP_TEXTURE_REPLACE);
	// 	glupPrivateTexCoord1d(0.0);
	// 	// glupSetMeshWidth(2.);

	// 	glupBegin(GLUP_LINES);
	// 	for(auto he : tet.iter_halfedges()) {
	// 		auto f = he.from().pos();
	// 		auto t = he.to().pos();

	// 		glupPrivateVertex3dv(GEO::vec3(f.x, f.y, f.z).data()); // draw first vertex
	// 		glupPrivateVertex3dv(GEO::vec3(t.x, t.y, t.z).data()); // draw second vertex
	// 	}
	// 	glupDisable(GLUP_TEXTURING);
	// 	glupEnd();
	// }

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

	ImGui::Checkbox("Display grid", &display_grid_);
	ImGui::Checkbox("Display axes", &display_axes_);
}

void App::GL_initialize() {
    SimpleMeshApplicationExt::GL_initialize();
    // init_rgba_colormap("labeling",6,1,labeling_colors_.as_chars());
    // init_rgba_colormap("validity",2,1,validity_colors_.as_chars());
    // state_transition(state_); // not all state_transition() code has been executed if GL was not initialized (in particular because missing colormaps)
}

// TODO move in gui_base !
// index_t App::pickup_cell_edge(GEO::vec3 p0, index_t c_idx) {
// 	// Search nearest edge
// 	double min_dist = std::numeric_limits<double>().max();
// 	index_t e_idx = -1;

// 	// Search nearest edge
// 	for (int i = 0; i < mesh_.cells.nb_edges(c_idx); i++) {
// 		// Get points from current edge
// 		index_t v0_idx = mesh_.cells.edge_vertex(c_idx, i, 0);
// 		index_t v1_idx = mesh_.cells.edge_vertex(c_idx, i, 1);
// 		GEO::vec3 &p1 = mesh_.vertices.point(v0_idx);
// 		GEO::vec3 &p2 = mesh_.vertices.point(v1_idx);

// 		// Compute dist from mouse point to edge points
// 		double dist = length(cross(p0 - p1, p0 - p2)) / length(p2 - p1);

// 		// Keep min dist
// 		if (dist < min_dist) {
// 			min_dist = dist;
// 			e_idx = i;
// 			posA = p1;
// 			posB = p2;
// 		}
// 	}

// 	return e_idx;
// }
index_t App::pickup_cell_edge(GEO::vec3 p0, index_t c_idx) {
	// Search nearest edge
	double min_dist = std::numeric_limits<double>().max();
	index_t e_idx = -1;

	// Search nearest edge
	for (int i = 0; i < mesh_.cells.nb_edges(c_idx); i++) {
		// Get points from current edge
		index_t v0_idx = mesh_.cells.edge_vertex(c_idx, i, 0);
		index_t v1_idx = mesh_.cells.edge_vertex(c_idx, i, 1);
		GEO::vec3 &p1 = mesh_.vertices.point(v0_idx);
		GEO::vec3 &p2 = mesh_.vertices.point(v1_idx);

		// Compute dist from mouse point to edge points
		double dist = distance(p0, (p1 + p2) / 2.);

		// Keep min dist
		if (dist < min_dist) {
			min_dist = dist;
			e_idx = i;
			posA = p1;
			posB = p2;
		}
	}

	return e_idx;
}

std::tuple<index_t, index_t> App::pickup_cell_facet(GEO::vec3 p0, index_t c_idx) {
	// Search if point is on facet
	double min_dist = std::numeric_limits<double>().max();
	index_t f_idx = NO_FACET;
	index_t lf_idx = NO_FACET;

	for (index_t lf = 0; lf < mesh_.cells.nb_facets(c_idx); lf++) {

		auto a = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 0));
		auto b = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 1));
		auto c = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 2));
		auto d = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 3));

		auto bary = (a + b + c + d) / 4.;

		auto n = normalize(cross(b - a, c - b));
		double dist = dot(p0 - a, n);
		
		if (std::abs(dot(normalize(p0 - a), normalize(bary - a))) < 1e-4)
			continue;

		if (dist < min_dist) {
			min_dist = dist;
			f_idx = mesh_.cells.facet(c_idx, lf);
			lf_idx = lf;
		}
	}

	return std::tuple<index_t, index_t>(f_idx, lf_idx);
}

void loop_cut(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Facet&)> f) {
	assert(hex.connected());

	std::vector<bool> visited(24 * hex.ncells(), false);
	std::queue<int> q;
	q.push(start_he);
	visited[start_he] = true;


	while (!q.empty()) {

		auto he_idx = q.front();
		auto he = Volume::Halfedge(hex, he_idx);

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

		// Propagate to front
		opp_c = he.opposite_f().next().next().opposite_f().opposite_c();

		if (opp_c.active()) {
			auto n_he = opp_c;

			if (!visited[n_he]) {
				q.push(n_he);
				visited[n_he] = true;
			}
		} 

		// Process current
		auto facet = he.opposite_f().facet();
		f(facet);
		
	}
}

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
	if (gui_mode == Camera || is_loading)
		return;

	// Try to pick cell
	hovered_cell = pick(MESH_CELLS);


	// If a cell is hovered, try to pick cell edge / cell facet 
	if (is_cell_hovered()) {
		index_t e_idx = pickup_cell_edge(picked_point_, hovered_cell);
		auto [f_idx, lf_idx] = pickup_cell_facet(picked_point_, hovered_cell);

		hovered_edge = e_idx;
		hovered_facet = f_idx;
		hovered_lfacet = lf_idx;
	}

	// TODO check MESH TYPE !!
	if (gui_mode == LoopCutPad) {

		hovered_path.clear();

		if (hex.connected() && is_cell_hovered()) {

			Volume::Cell um_c(hex, hovered_cell);
			// auto he = um_c.halfedge(he_n);
			// posAb = he.from().pos();
			// posBb = he.to().pos();
			// Quad3 q = he.facet();
			// posN = q.normal();


			
			Volume::Halfedge hovered_he(hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(hovered_edge, hovered_lfacet)));
			
			if (hovered_he.active()) {
				loop_cut(hex, hovered_he, [&](Volume::Halfedge &he) {
					UM::vec3 a = he.from().pos();
					UM::vec3 b = he.to().pos();
					hovered_path.push_back(a);
					hovered_path.push_back(b);
				});
			}
		}
	}
}

void App::mouse_button_callback(int button, int action, int mods, int source) {

	if (gui_mode == Camera) {
    	SimpleMeshApplication::mouse_button_callback(button,action,mods,source);
		return;
	}

	// If left click
    if (gui_mode == Painting && action == EVENT_ACTION_DOWN && button == 0) {

		index_t f_idx = pick(MESH_FACETS);
		if (f_idx != NO_FACET && f_idx < mesh_.facets.nb()) {
			std::cout << "facet: " << f_idx << std::endl;
			GEO::Attribute<GEO::Numeric::uint32> flag(
				mesh_.facets.attributes(), "flag"
			);

			flag[f_idx] = paint_value;

			index_t f_idx = pick(MESH_FACETS);
			std::cout << "pos: " << picked_point_ << std::endl;
			click_pos = picked_point_;

		}

    }
	else if (gui_mode == Painting && action == EVENT_ACTION_UP && button == 0) {
		// um_bindings::um_tet_from_geo_mesh(mesh_, tet);
	}
	else if (gui_mode == LoopCutPad && action == EVENT_ACTION_DOWN && button == 0) {

		// TODO move because it will be useful in other mode than Layer
		selected_vertex = hovered_vertex;
		selected_edge = hovered_edge;
		selected_facet = hovered_facet;
		selected_lfacet = hovered_lfacet;
		selected_cell = hovered_cell;

		selected_path = hovered_path;

		// Test extract layer
		// Quads q_out;
		// Volume::Cell um_c(hex, hovered_cell);
		// Volume::Halfedge start_he(hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(selected_edge, selected_lfacet)));
		// loop_cut(hex, start_he, [&](UM::Volume::Facet &f) {
			
		// 	int v_off = q_out.points.create_points(4);
		// 	q_out.points[v_off] = f.vertex(0).pos();
		// 	q_out.points[v_off + 1] = f.vertex(1).pos();
		// 	q_out.points[v_off + 2] = f.vertex(2).pos();
		// 	q_out.points[v_off + 3] = f.vertex(3).pos();

		// 	int f_off = q_out.create_facets(1);
		// 	q_out.vert(f_off, 0) = v_off;
		// 	q_out.vert(f_off, 1) = v_off + 1;
		// 	q_out.vert(f_off, 2) = v_off + 2;
		// 	q_out.vert(f_off, 3) = v_off + 3;

		// });
		// write_by_extension("result.geogram", q_out);
	}



}

void App::key_callback(int key, int scancode, int action, int mods) {
	
	std::cout << "key pressed: " << key << std::endl;

	// Ctrl
	if (action == EVENT_ACTION_DOWN && key == 341) {
		switch_mode = gui_mode;
		gui_mode = Camera;
	} else if (action == EVENT_ACTION_UP && key == 341) {
		gui_mode = switch_mode;
		switch_mode = Camera;
	}

	// Validate layer
	if (gui_mode == LoopCutPad && (key == 257 || key == 335) && action == EVENT_ACTION_DOWN && is_cell_selected()) {
		
		CellFacetAttribute<bool> pad_face(hex);

		Volume::Cell um_c(hex, selected_cell);
		Volume::Halfedge start_he(hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(selected_edge, selected_lfacet)));

		loop_cut(hex, start_he, [&](UM::Volume::Facet &f) {
			pad_face[f] = true;
		});

		BenjaminAPI::pad(hex, pad_face);

		// TODO here use hexboundary
		um_bindings::geo_mesh_from_um_hex(hex, mesh_);
		mesh_gfx_.set_mesh(&mesh_);

		selected_path.clear();
		write_by_extension("padded.geogram", hex);
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

bool App::load(const std::string& filename) {
    
	mesh_gfx_.set_mesh(nullptr);
    mesh_.clear(false, true);

    MeshIOFlags flags;
	
	is_loading = true;

    if(!mesh_load(filename, mesh_, flags)) {
		// TODO clear all here too !
        return false;
    }

	// Reset all values
	// gui_mode = Camera;
	// Clear selections
	reset_hovered_selected();
	// Clear path
	hovered_path.clear();
	selected_path.clear();

	// Clear UM meshes
	tet.clear();
	hex.clear();

	is_loading = false;

	// Init UM tet from GEO mesh
	um_bindings::um_tet_from_geo_mesh(mesh_, tet);
	tet.connect();



	// Display info
    mesh_gfx_.set_animate(false);
    mesh_.vertices.set_dimension(3);
    mesh_.vertices.set_double_precision(); // added
    mesh_gfx_.set_mesh(&mesh_);
    current_file_ = filename;

	// TODO see if useless
    labeling_visu_mode_transition();

    clear_scene_overlay();

	std::cout << "Mesh loaded !" << std::endl;

    return true;
}


void App::draw_object_properties() {

	ImGui::Text("Hovered vertex: %i", hovered_vertex);
	ImGui::Text("Hovered cell edge: %i", hovered_edge);
	ImGui::Text("Hovered cell facet: %i", hovered_facet);
	ImGui::Text("Hovered cell local facet: %i", hovered_lfacet);
	ImGui::Text("Hovered cell: %i", hovered_cell);
	// ImGui::Text("Hen: %i", he_n);

	ImGui::Checkbox("Show vertices", &show_vertices_);
	ImGui::Checkbox("Show surface", &show_surface_);
	ImGui::Checkbox("Show volume", &show_volume_);

	// if(ImGui::Button("Next!")) {
	// 	he_n++;
	// }

	if(ImGui::Button("Camera")) {
		gui_mode = Camera;
	}

	if(ImGui::Button("Hover")) {
		gui_mode = Hover;
	}

	ImGui::Separator();

	if(ImGui::Button("Paint X")) {
		paint_value = 0;
		gui_mode = Painting;
	}
	if(ImGui::Button("Paint Y")) {
		paint_value = 1;
		gui_mode = Painting;
	}
	if(ImGui::Button("Paint Z")) {
		paint_value = 2;
		gui_mode = Painting;
	}	
	if(ImGui::Button("Paint -X")) {
		paint_value = 3;
		gui_mode = Painting;
	}
	if(ImGui::Button("Paint -Y")) {
		paint_value = 4;
		gui_mode = Painting;
	}
	if(ImGui::Button("Paint -Z")) {
		paint_value = 5;
		gui_mode = Painting;
	}	

	if(ImGui::Button("Layer view !")) {
		gui_mode = LoopCutPad;
	}	

	ImGui::Separator();

	// TODO Add checking of mesh type ! Maybe a class that derives from TetBoundary to rewrite with check() override...
	// TODO Add view mode for mesh type to switch view

	if (ImGui::Button("Compute flags !")) {
		
		// Compute flag on tet and tri
		TetBoundary tet_bound(tet);
		
		// To GEO mesh
		um_bindings::geo_mesh_from_tetboundary(tet_bound, mesh_);

		UM::CellFacetAttribute<int> tet_flag(tet, -1);
		UM::FacetAttribute<int> tri_flag(tet_bound.tri, -1);

		// Compute flag
		algo::naive_tag(tet, tet_flag);
		// Transfert flag from tet to tri for display
		tet_bound.set_attribute_to_surface(tet_flag, tri_flag);
		// Update GEO mesh attribute "flag"
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_CELL_FACETS>(tet, tet_flag.ptr, "tet_flag", mesh_);
		um_bindings::geo_attr_from_um_attr2<GEO::MESH_FACETS>(tet_bound.tri, tri_flag.ptr, "flag", mesh_);

		labeling_visu_mode_transition();
		show_surface_ = true;
		show_volume_ = false;
	}

	ImGui::Separator();
	int nhex_wanted = 3000;

	if (ImGui::Button("Polycubify !")) {

		// Get UM cell facet attribute tet_flag from GEO mesh
		UM::CellFacetAttribute<int> tet_flag(tet, -1);
		um_bindings::um_attr_from_geo_attr<GEO::MESH_CELL_FACETS>(mesh_, "tet_flag", tet, tet_flag.ptr);

		// UM::Hexahedra hex;

		BenjaminAPI::polycubify(tet, tet_flag, hex, nhex_wanted);
		std::cout << "polycubify..." << std::endl;

		HexBoundary hex_bound(hex);
		// Replace current GEO mesh by UM Hex
		um_bindings::geo_mesh_from_hexboundary(hex_bound, mesh_);
		mesh_gfx_.set_mesh(&mesh_);

		labeling_visu_mode_transition();
		show_surface_ = false;
		// show_surface_borders_ = true;
		show_volume_ = true;
		// show_hexes_ = true;

		// GEO::Attribute<int> geo_attr(
		// 	mesh_.cells.attributes(), "test"
		// );

		// for (int i = 0; i < geo_attr.size(); i++)
		// 	geo_attr[i] = rand() % 100;

		// // attribute_ = "cells.test";
		// // attribute_subelements_ = MESH_CELLS;
		// // attribute_name_ = "test";
		// // attribute_min_ = 0;
		// // attribute_max_ = 100;

		// // attribute_ = "cell_facets.tet_flag";
		// // attribute_subelements_ = MESH_CELL_FACETS;
		// // attribute_name_ = "tet_flag";
		// // attribute_min_ = -1;
		// // attribute_max_ = 5;

		write_by_extension("input.geogram", tet);
		write_by_extension("polycubify_hex.geogram", hex);
		write_by_extension("polycubify_quad.geogram", hex_bound.quad);
		mesh_save(mesh_, "polycubify2.geogram");

	}

}



// void App::draw_object_properties() {

// 	ImGui::Checkbox("Show vertices", &show_vertices_);
// 	ImGui::Checkbox("Show surface", &show_surface_);
// 	ImGui::Checkbox("Show volume", &show_volume_);

// 	if(ImGui::Button("Camera tool")) {
// 		gui_mode = Camera;
// 	}

// 	ImGui::Separator();

// 	if(ImGui::Button("Paint X")) {
// 		paint_value = 0;
// 		gui_mode = Painting;
// 	}
// 	if(ImGui::Button("Paint Y")) {
// 		paint_value = 1;
// 		gui_mode = Painting;
// 	}
// 	if(ImGui::Button("Paint Z")) {
// 		paint_value = 2;
// 		gui_mode = Painting;
// 	}	
// 	if(ImGui::Button("Paint -X")) {
// 		paint_value = 3;
// 		gui_mode = Painting;
// 	}
// 	if(ImGui::Button("Paint -Y")) {
// 		paint_value = 4;
// 		gui_mode = Painting;
// 	}
// 	if(ImGui::Button("Paint -Z")) {
// 		paint_value = 5;
// 		gui_mode = Painting;
// 	}	

// 	ImGui::Separator();

// 	// TODO Add checking of mesh type ! Maybe a class that derives from TetBoundary to rewrite with check() override...
// 	// TODO Add view mode for mesh type to switch view

// 	if(ImGui::Button("Compute flags !")) {
		
// 		um_bindings::combinatorial_update(mesh_, [&](TetBoundary& tet_bound, SurfaceAttributes& surf_attr, VolumeAttributes& vol_attr) {

// 			// Declare flag attribute
// 			UM::CellFacetAttribute<int> tet_flag(tet_bound.tet, -1);
// 			UM::FacetAttribute<int> tri_flag(tet_bound.tri, -1);

// 			// Compute flag on tet and tri
// 			algo::naive_tag(tet_bound.tet, tet_flag);
// 			// Transfert cell facet attribute of tet to tri surface
// 			tet_bound.set_attribute_to_surface(tet_flag, tri_flag);

// 			// Return attributes i would like to transfert to mesh
// 			return std::pair<SurfaceAttributes, VolumeAttributes>{{{}, {{"flag", tri_flag.ptr}}, {}}, {{}, {}, {{"flag", tet_flag.ptr}}, {}}};
// 		});
		
// 		labeling_visu_mode_transition();
// 	}



// 	if(ImGui::Button("Integrate !")) {

// 		um_bindings::combinatorial_update(mesh_, [&](TetBoundary& tet_bound, SurfaceAttributes& surf_attr, VolumeAttributes& vol_attr) {

// 			write_by_extension("integrated_tri.geogram", tet_bound.tri);

// 			UM::PointAttribute<UM::vec3> U(tet_bound.tet);

// 			UM::FacetAttribute<int> tri_flag("flag", surf_attr, tet_bound.tri);
// 			UM::CellFacetAttribute<int> tet_flag("flag", vol_attr, tet_bound.tet);


// 			BenjaminAPI::integrate(tet_bound.tet, tet_flag, U);
// 			// Update vertices pos
// 			for (auto v : tet_bound.tet.iter_vertices()) {
// 				v.pos() = U[v];
// 			}
// 			for (auto v : tet_bound.tri.iter_vertices()) {
// 				v.pos() = U[v];
// 			}

// 			// write_by_extension("integrated_tet.geogram", tet_bound.tet, {{}, {}, {{"flag", tet_flag.ptr}}, {}});
// 			// write_by_extension("integrated_tri.geogram", tet_bound.tri, {{}, {{"flag", tri_flag.ptr}}, {}});


// 			// Transfert cell facet attribute of tet to tri surface
// 			// UM::FacetAttribute<int> tri_flag(tet_bound.tri, -1);
// 			// tet_bound.set_attribute_to_surface(tet_flag, tri_flag);

// 			// write_by_extension("integrated_tet.geogram", tet_bound.tet);
// 			// write_by_extension("integrated_tri.geogram", tet_bound.tri);

// 			std::cout << "integrated !!!" << std::endl;

// 			// TODO should return tri_flag for refreshing view ! 
// 			// return std::pair<SurfaceAttributes, VolumeAttributes>{{{}, {}, {}}, {{}, {}, {}, {}}};
// 			return std::pair<SurfaceAttributes, VolumeAttributes>{{{}, {{"flag", tri_flag.ptr}}, {}}, {{}, {}, {{"flag", tet_flag.ptr}}, {}}};
// 		});

// 		labeling_visu_mode_transition();
// 	}


// 	// if(ImGui::Button("Randomize !")) {

// 	// 	um_bindings::combinatorial_update(mesh_, [&](TetBoundary& tet_bound, SurfaceAttributes& surf_attr, VolumeAttributes& vol_attr) {

// 	// 		for (auto v : tet_bound.tri.iter_vertices()) {
// 	// 			v.pos() += {rand()%100/10000., rand()%100/10000., rand()%100/10000.};
// 	// 		}

// 	// 		UM::CellFacetAttribute<int> tet_flag(tet_bound.tet, -1);
// 	// 		UM::FacetAttribute<int> tri_flag(tet_bound.tri, -1);

// 	// 		return std::pair<SurfaceAttributes, VolumeAttributes>{{{}, {{"flag", tri_flag.ptr}}, {}}, {{}, {}, {{"flag", tet_flag.ptr}}, {}}};
// 	// 	});
// 	// }

// 	// if(ImGui::Button("Compute flags !")) {

// 	// 	// Extract surface
// 	// 	TetBoundary tet_b(tet);

// 	// 	// Compute flag on tet and tri
// 	// 	UM::CellFacetAttribute<int> tet_flag(tet_b.tet, -1);
// 	// 	algo::naive_tag(tet, tet_flag);		

// 	// 	// Transfert cell facet attribute of tet to tri surface
// 	// 	UM::FacetAttribute<int> tri_flag(tet_b.tri, -1);
// 	// 	tet_b.set_attribute_to_surface(tet_flag, tri_flag);

// 	// 	// UM surface to mesh
// 	// 	// um_bindings::geo_mesh_from_um_tri(tet_b.tri, mesh_);
// 	// 	// um_bindings::geo_attr_from_um_attr(tet_b.tri, "flag", tri_flag, mesh_);

		
// 	// 	UM::PointAttribute<UM::vec3> U(tet);
// 	// 	BenjaminAPI::integrate(tet, tet_flag, U);
// 	// 	for (auto v : tet.iter_vertices()) {
// 	// 		v.pos() = U[v];
// 	// 	}
// 	// 	tet_b.set(tet);
// 	// 	um_bindings::geo_mesh_from_um_tri(tet_b.tri, mesh_);
// 	// 	um_bindings::geo_attr_from_um_attr(tet_b.tri, "flag", tri_flag, mesh_);

// 	// 	// std::cout << "integrated..." << std::endl;
// 	// 	// UM::write_by_extension("integrated.geogram", tet, {});
// 	// 	// UM::write_by_extension("integrated_tri.geogram", tet_b.tri, {{}, {{"flag", tri_flag.ptr}}, {}});

// 	// 	// labeling_visu_mode_transition();

// 	// 	/*
// 	// 	UM::Hexahedra hex;
// 	// 	BenjaminAPI::polycubify(tet, tet_flag, hex, 1000);
// 	// 	std::cout << "polycubify..." << std::endl;

// 	// 	UM::write_by_extension("polycubify.geogram", hex, {});
// 	// 	*/

// 	// 	// std::cout << "Mesh writed by ultimaille" << std::endl;
// 	// }

// 	// if(ImGui::Button("Polycubify !")) {

// 	// 	std::cout << "polycubify..." << std::endl;

// 	// 	// Need to retransform mesh to tet

// 	// 	UM::CellFacetAttribute<int> tet_flag(tet);
// 	// 	um_bindings::um_attr_from_geo_attr(mesh_, "flag", tet, tet_flag);
// 	// 	std::cout << "polycubify..." << std::endl;
		
// 	// 	for (auto v : tet.iter_vertices()) {
// 	// 		std::cout << "vertices: " << v << std::endl;
// 	// 	}

// 	// 	UM::Hexahedra hex;
// 	// 	BenjaminAPI::polycubify(tet, tet_flag, hex, 1000);
// 	// 	std::cout << "polycubify..." << std::endl;

// 	// 	UM::write_by_extension("polycubify.geogram", hex, {});
// 	// 	std::cout << "polycubified !" << std::endl;
// 	// }

// }

std::string App::supported_write_file_extensions() {
    return SimpleMeshApplication::supported_write_file_extensions() + ";txt"; // add .txt in supported write file extensions
}

void App::labeling_visu_mode_transition() {
    if(colormaps_.empty()) {
        // GL is not initialized yet
        // the state transition will be triggered later, in GL_initialize()
        std::cout << "empty colormap ! " << std::endl;
        return;
    }

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