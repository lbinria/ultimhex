#include "path_constraint_padding_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"


void PathConstraintPaddingTool::puff_view(Volume::Halfedge &h) {
	helpers::get_halfedge_layers(ctx.hex_bound->hex, layer);

	CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
	helpers::puff(ctx.hex_bound->hex, h, layer, selected);

	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	// ctx.recompute_hex(selected);
	ctx.view.switch_to_surface_select_mode();
	ctx.gui_mode = PathConstraintPadding;
}

void PathConstraintPaddingTool::reset_attr() {
	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);


	// for (int i = 0; i < v_cell_facets_hovered_attr.size(); i++) {
	// 	cell_facets_hovered_attr[i] = v_cell_facets_hovered_attr[i];
	// 	auto qf = ctx.hex_bound->quad_facet(i);
	// 	if (qf >= 0 && qf < ctx.mesh_.facets.nb())
	// 		hovered_attr[qf] = v_hovered_attr[qf];
	// }

	// TODO reset ?

	for (int f = 0; f < v_cell_facets_hovered_attr.size(); f++) {
		cell_facets_hovered_attr[f] = v_cell_facets_hovered_attr[f];
		// auto qf = ctx.hex_bound->quad_facet(f);
		// if (qf >= 0 && qf < ctx.mesh_.facets.nb())
		// 	hovered_attr[qf] = 2;
	}
}

bool PathConstraintPaddingTool::draw_object_properties() {

	auto start_h = *ctx.hex_bound->hex.iter_halfedges().begin();

	

	if (ImGui::Button("Puff pastry 0##btn_path_constraint_padding_tool_puff_0", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		puff_view(start_h);
		reset_attr();
		return true;
	}

	if (ImGui::Button("Puff pastry 1##btn_path_constraint_padding_tool_puff_1", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		auto h = start_h.opposite_f().prev();
		puff_view(h);
		reset_attr();
		return true;
	}

	if (ImGui::Button("Puff pastry 2##btn_path_constraint_padding_tool_puff_2", ImVec2(ImGui::GetContentRegionAvail().x, 0))) {
		auto h = start_h.opposite_f().prev().opposite_f().prev();
		puff_view(h);
		reset_attr();
		return true;
	}

	if (ImGui::SliderFloat3("Plane center##path_constraint_padding_tool_plane_center", v_plane_center, -1., 1.)) {
		plane_center = GEO::vec3(v_plane_center[0], v_plane_center[1], v_plane_center[2]);
		ctx.gui_mode = PathConstraintPadding;
	}
	if (ImGui::SliderAngle("Plane center##path_constraint_padding_tool_plane_rotx", &v_plane_rot_x, 0.f, 180.f)) {
		plane_rot = GEO::vec3(v_plane_rot_x, v_plane_rot_y, v_plane_rot_z);
		ctx.gui_mode = PathConstraintPadding;
	}
	if (ImGui::SliderAngle("Plane center##path_constraint_padding_tool_plane_roty", &v_plane_rot_y, 0.f, 180.f)) {
		plane_rot = GEO::vec3(v_plane_rot_x, v_plane_rot_y, v_plane_rot_z);
		ctx.gui_mode = PathConstraintPadding;
	}
	if (ImGui::SliderAngle("Plane center##path_constraint_padding_tool_plane_rotz", &v_plane_rot_z, 0.f, 180.f)) {
		plane_rot = GEO::vec3(v_plane_rot_x, v_plane_rot_y, v_plane_rot_z);
		ctx.gui_mode = PathConstraintPadding;
	}

	if (ImGui::Button("Slice##path_constraint_padding_tool_slice")) {
		

		slice();
		return true;
	}

	return false;
}

void PathConstraintPaddingTool::slice() {
		// Search for cell that contains the point
		Volume::Cell found_c(ctx.hex_bound->hex, -1);

		for (auto c : ctx.hex_bound->hex.iter_cells()) {
			double volume = 0;
			for (auto f : c.iter_facets()) {
				Pyramid p {
					f.vertex(0).pos(),
					f.vertex(1).pos(),
					f.vertex(2).pos(),
					f.vertex(3).pos(),
					um_bindings::um_vec(plane_center)
				};
				volume += abs(p.volume());
			}
			// volume *= .5;
			std::cout << "Volume from point: " << volume << std::endl;
			std::cout << "Hex volume: " << c.geom<Hexahedron>().volume() << std::endl;

			if (abs(volume - c.geom<Hexahedron>().volume()) < 1e-8) {
				found_c = c;
				break;
			}
		}

		if (!found_c.active()) {
			std::cout << "No cell found" << std::endl;
			return;
		}

		// Compute plane normal
		Quad3 pq { 
			um_bindings::um_vec(v_plane_corners[0]),
			um_bindings::um_vec(v_plane_corners[1]),
			um_bindings::um_vec(v_plane_corners[2]),
			um_bindings::um_vec(v_plane_corners[3]),
		};

		std::cout << "Plane normal: " << pq.normal() << std::endl;

		Volume::Facet found_f(ctx.hex_bound->hex, -1);
		double min_angle = std::numeric_limits<double>::max();

		for (auto f : found_c.iter_facets()) {
			Quad3 fq = f.geom<Quad3>();
			std::cout << "Facet normal: " << fq.normal() << std::endl;
			double angle = 1. - fq.normal() * pq.normal();
			std::cout << "Angle: " << angle << std::endl;

			if (angle < min_angle) {	
				min_angle = angle;
				found_f = f;
			}
		}

		// Get layer from facet
		helpers::get_halfedge_layers(ctx.hex_bound->hex, layer);
		int l = layer[found_f.halfedge(0).opposite_f().next()];


		// Get arbitrary facet
		GEO::Attribute<int> cell_facets_hovered_attr(
			ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
		);

		for (auto f : ctx.hex_bound->hex.iter_facets()) {
			if (layer[f.halfedge(0).opposite_f().next()] == l) {
				cell_facets_hovered_attr[f] = 2;
			}
		}


		// Update view
		CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
		for (auto f : ctx.hex_bound->hex.iter_facets()) {
			if (cell_facets_hovered_attr[f] >= 2) {
				selected[f] = true;
			}
		}

		// Extract surface from selected facets
		ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		ctx.view.switch_to_surface_select_mode();


}

void PathConstraintPaddingTool::draw_viewer_properties() {

}

void PathConstraintPaddingTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
	glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec4f(1., 0., 0., 1.).data());
	glupSetColor4fv(GLUP_BACK_COLOR, GEO::vec4f(1., 0., 0., 1.).data());

	glupBegin(GLUP_POINTS);
	glupPrivateVertex3dv(plane_center.data());

	GEO::vec3 corners[] = {
		GEO::vec3(-0.5, -0.5, 0.),
		GEO::vec3(0.5, -0.5, 0.),
		GEO::vec3(0.5, 0.5, 0.),
		GEO::vec3(-0.5, 0.5, 0.)
	};

	for (int lc = 0; lc < 4; lc++) {
		double rad = plane_rot.z;

		GEO::vec3 r = corners[lc];
		
		GEO::vec3 r0;
		r0.x = r.x * cos(rad) - r.y * sin(rad);
		r0.y = r.x * sin(rad) + r.y * cos(rad);
		r0.z = r.z;
		r = r0;

		rad = plane_rot.y;
		GEO::vec3 r1;
		r1.x = r.x * cos(rad) + r.z * sin(rad);
		r1.y = r.y;
		r1.z = -r.x * sin(rad) + r.z * cos(rad);
		r = r1;

		rad = plane_rot.x;
		
		GEO::vec3 r2;
		r2.x = r.x;
		r2.y = r.y * cos(rad) - r.z * sin(rad);
		r2.z = r.y * sin(rad) + r.z * cos(rad);
		r = r2;

		v_plane_corners[lc] = r;
		glupPrivateVertex3dv((plane_center + r).data());
	}

	glupEnd();

}

void PathConstraintPaddingTool::hover_callback(double x, double y, int source) {
	// TODO BIG COPY/PASTE FROM BLOC PAD TOOL
	// Refactor this
	if (/*step > 0 ||*/ !ctx.is_facet_hovered())
		return;


	std::vector<int> hovered_facets;

	// Chart selection mode
	// if (select_mode == 1) {
	// 	for (auto f : ctx.hex_bound->quad.iter_facets()) {
	// 		if (patches[ctx.hex_bound->hex_facet(f)] == patches[ctx.hex_bound->hex_facet(ctx.hovered_facet)])
	// 			hovered_facets.push_back(f);
	// 	}
	// } else {
		hovered_facets.push_back(ctx.hovered_facet);
	// }

	// Attribute hovered / selected, enable visualizing hovered / selected facets
	GEO::Attribute<int> hovered_attr(
		ctx.mesh_.facets.attributes(), "hovered"
	);
	GEO::Attribute<int> cell_facets_hovered_attr(
		ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	);

	// Reset last hovered facets as not hovered, if not selected
	for (auto last_hovered_f : last_hovered_facets) {
		if (last_hovered_f >= 0 && hovered_attr[last_hovered_f] == 1) {
			hovered_attr[last_hovered_f] = 0;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(last_hovered_f)] = 0;
		}
	}

	// Last hovered facet become current hovered facet
	last_hovered_facets = hovered_facets;

	for (auto f : hovered_facets) {
		// If not selected, facet is hovered
		if (hovered_attr[f] < 2) {
			hovered_attr[f] = 1;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 1;
		}

		if (ctx.left_mouse_pressed) {
			hovered_attr[f] = 2;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 2;
		} else if (ctx.right_mouse_pressed) {
			hovered_attr[f] = 0;
			cell_facets_hovered_attr[ctx.hex_bound->hex_facet(f)] = 0;
		}
	}

	v_cell_facets_hovered_attr = cell_facets_hovered_attr.get_vector();
}

void PathConstraintPaddingTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void PathConstraintPaddingTool::scroll_callback(double xoffset, double yoffset) {}

void PathConstraintPaddingTool::validate_callback() {

	if (mode == Select) {

		GEO::Attribute<int> cell_facets_hovered_attr(
			ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
		);
		CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
		for (auto f : ctx.hex_bound->hex.iter_facets()) {
			if (cell_facets_hovered_attr[f] >= 2) {
				selected[f] = true;
			}
		}

		// Extract surface from selected facets
		ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		ctx.view.switch_to_surface_select_mode();


		// GEO::Attribute<int> cell_facets_hovered_attr(
		// 	ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
		// );

		// // Extends the selected facets to their layer 
		// CellFacetAttribute<bool> selected(ctx.hex_bound->hex, false);
		// std::vector<int> selector;

		// for (auto h : ctx.hex_bound->hex.iter_halfedges()) {
		// 	// Test if facet is selected
		// 	auto hf = h.prev().opposite_f();
		// 	if (cell_facets_hovered_attr[hf.facet()] < 2)
		// 		continue;

		// 	for (auto inner_h : ctx.hex_bound->hex.iter_halfedges()) {
		// 		// Test if selected facet is in the same layer as current facet
		// 		auto inner_hf = inner_h.prev().opposite_f();
		// 		if (layer[inner_h] != layer[h])
		// 			continue;

		// 		// Add current facet to selected facets
		// 		selected[inner_hf.facet()] = true;
		// 		selector.push_back(inner_hf.facet());
		// 	}
		// }

		// // Intersect layers
		// // for (auto f : ctx.hex_bound->hex.iter_facets()) {
		// // 	if (!selected[f])
		// // 		continue;
		// // }
		// DisjointSet ds(ctx.hex_bound->hex.nfacets());
		// for (auto fi : selector) {
		// 	Volume::Facet f(ctx.hex_bound->hex, fi);

		// 	for (auto h : f.iter_halfedges()) {
		// 		auto ortho_h = h.opposite_f().next();

		// 		auto opp = h.opposite_f().opposite_c(); 
		// 		if (opp.active()) {
		// 			auto opp_h = opp.opposite_f().next().next();
		// 			auto opp_ortho_h = opp_h.opposite_f().next();

		// 			// Check if the facets are in the same layer
		// 			if (layer[ortho_h] != layer[opp_ortho_h])
		// 				continue;
					
		// 			// Check there is no intersection between two layers
		// 			bool intersect_found = false;
		// 			for (auto ah : h.iter_CCW_around_edge()) {
		// 				// If layer wasn't selected exit
		// 				if (!selected[ah.facet()])
		// 					continue;

		// 				// If orthogonal layer is different, there is an intersection
		// 				auto ortho_ah = ah.opposite_f().next();
		// 				if (layer[ortho_h] != layer[ortho_ah]) {
		// 					intersect_found = true;
		// 					break;
		// 				}
		// 			}
					
		// 			// If there is no intersection, merge the two facets
		// 			if (!intersect_found)
		// 				ds.merge(f, opp_h.facet());
		// 		}
		// 	}

		// }

		// std::vector<int> facet_groups;
		// int nsets = ds.get_sets_id(facet_groups);
		

		// {
		// Quads q;
		// FacetAttribute<int> bob(q, -1);
		// for (auto f : ctx.hex_bound->hex.iter_facets()) {
		// 	if (selected[f]) {
				
		// 		int foff = q.create_facets(1);
				
		// 		for (int lv = 0; lv < 4; lv++) {
		// 			int voff = q.points.create_points(1);
		// 			q.points[voff] = f.vertex(lv).pos();
		// 			q.vert(foff, lv) = voff;
		// 		}

		// 		bob[foff] = facet_groups[f];
		// 	}
		// }

		// write_by_extension("bobi.geogram", q, {{}, {{"group", bob.ptr}}, {}});
		// }
		// std::set<int> groups;
		// for (auto f : selector) {
		// 	if (cell_facets_hovered_attr[f] >= 2)
		// 		groups.insert(facet_groups[f]);
		// }

		// CellFacetAttribute<bool> pad_facets(ctx.hex_bound->hex, false);
		// // Keep only sides where facets was selected 
		// {
		// Quads q;

		// for (auto f : ctx.hex_bound->hex.iter_facets()) {
		// 	if (groups.find(facet_groups[f]) == groups.end())
		// 		continue;

		// 	int foff = q.create_facets(1);
			
		// 	for (int lv = 0; lv < 4; lv++) {
		// 		int voff = q.points.create_points(1);
		// 		q.points[voff] = f.vertex(lv).pos();
		// 		q.vert(foff, lv) = voff;
		// 	}
		// 	pad_facets[f] = true;	
		// }
		// write_by_extension("layer_intersections.geogram", q, {{}, {}, {}});
		// }

		

		// // Extract surface from selected facets
		// // ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, selected);
		// ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex, pad_facets);
		// um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);
		// ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
		// ctx.view.switch_to_surface_select_mode();

		mode = Preview;
	}
	if (mode == Preview) {
		
	}
}

void PathConstraintPaddingTool::escape_callback() {
	if (mode == Preview) {
		clear();
	}
}

bool PathConstraintPaddingTool::is_compatible() { 
	return true;
}