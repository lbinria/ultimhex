#include "path_constraint_padding_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"
#include "../helpers.h"


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


		Volume::Facet found_f(ctx.hex_bound->hex, -1);
		double min_angle = std::numeric_limits<double>::max();

		for (auto f : found_c.iter_facets()) {
			Quad3 fq = f.geom<Quad3>();
			double angle = 1. - fq.normal() * pq.normal();

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
		ctx.view.switch_to_surface_select_mode();

		ctx.view.show_volume_ = true;
		ctx.view.cells_shrink_ = 0.4f;
		ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

}

void PathConstraintPaddingTool::reset_attr() {
	// // Attribute hovered / selected, enable visualizing hovered / selected facets
	// GEO::Attribute<int> hovered_attr(
	// 	ctx.mesh_.facets.attributes(), "hovered"
	// );
	// GEO::Attribute<int> cell_facets_hovered_attr(
	// 	ctx.mesh_.cell_facets.attributes(), "cell_facets_hovered"
	// );


	// // for (int i = 0; i < v_cell_facets_hovered_attr.size(); i++) {
	// // 	cell_facets_hovered_attr[i] = v_cell_facets_hovered_attr[i];
	// // 	auto qf = ctx.hex_bound->quad_facet(i);
	// // 	if (qf >= 0 && qf < ctx.mesh_.facets.nb())
	// // 		hovered_attr[qf] = v_hovered_attr[qf];
	// // }

	// // TODO reset ?

	// for (int f = 0; f < v_cell_facets_hovered_attr.size(); f++) {
	// 	cell_facets_hovered_attr[f] = v_cell_facets_hovered_attr[f];
	// 	// auto qf = ctx.hex_bound->quad_facet(f);
	// 	// if (qf >= 0 && qf < ctx.mesh_.facets.nb())
	// 	// 	hovered_attr[qf] = 2;
	// }
}

bool PathConstraintPaddingTool::draw_object_properties() {

	// auto start_h = *ctx.hex_bound->hex.iter_halfedges().begin();

	if (ImGui::SliderFloat3("Plane center##path_constraint_padding_tool_plane_center", v_plane_center, -1., 1.)) {
		plane_center = GEO::vec3(v_plane_center[0], v_plane_center[1], v_plane_center[2]);
		ctx.gui_mode = PathConstraintPadding;
	}
	if (ImGui::SliderAngle("Plane center##path_constraint_padding_tool_plane_rotx", &v_plane_rot_x, 0.f, 360.f)) {
		plane_rot = GEO::vec3(v_plane_rot_x, v_plane_rot_y, v_plane_rot_z);
		ctx.gui_mode = PathConstraintPadding;
	}
	if (ImGui::SliderAngle("Plane center##path_constraint_padding_tool_plane_roty", &v_plane_rot_y, 0.f, 360.f)) {
		plane_rot = GEO::vec3(v_plane_rot_x, v_plane_rot_y, v_plane_rot_z);
		ctx.gui_mode = PathConstraintPadding;
	}
	if (ImGui::SliderAngle("Plane center##path_constraint_padding_tool_plane_rotz", &v_plane_rot_z, 0.f, 360.f)) {
		plane_rot = GEO::vec3(v_plane_rot_x, v_plane_rot_y, v_plane_rot_z);
		ctx.gui_mode = PathConstraintPadding;
	}

	if (ImGui::Button("Slice##path_constraint_padding_tool_slice")) {
		slice();
		return true;
	}

	return false;
}

void PathConstraintPaddingTool::draw_viewer_properties() {

}

static GEO::vec3 rotate(GEO::vec3 v, GEO::vec3 rot) {
	double rad = rot.z;
	GEO::vec3 r0;
	r0.x = v.x * cos(rad) - v.y * sin(rad);
	r0.y = v.x * sin(rad) + v.y * cos(rad);
	r0.z = v.z;

	v = r0;

	rad = rot.y;
	GEO::vec3 r1;
	r1.x = v.x * cos(rad) + v.z * sin(rad);
	r1.y = v.y;
	r1.z = -v.x * sin(rad) + v.z * cos(rad);
	v = r1;

	rad = rot.x;
	
	GEO::vec3 r2;
	r2.x = v.x;
	r2.y = v.y * cos(rad) - v.z * sin(rad);
	r2.z = v.y * sin(rad) + v.z * cos(rad);
	v = r2;

	return v;
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
		v_plane_corners[lc] = rotate(corners[lc], plane_rot);
		glupPrivateVertex3dv((plane_center + v_plane_corners[lc]).data());
	}

	glupEnd();

	glupBegin(GLUP_LINES);
	for (int lc = 0; lc < 4; lc++) {
		glupPrivateVertex3dv((plane_center + v_plane_corners[lc]).data());
		glupPrivateVertex3dv((plane_center + v_plane_corners[(lc + 1) % 4]).data());
	}
	glupEnd();

	Quad3 pq {
		um_bindings::um_vec(v_plane_corners[0]),
		um_bindings::um_vec(v_plane_corners[1]),
		um_bindings::um_vec(v_plane_corners[2]),
		um_bindings::um_vec(v_plane_corners[3])
	};

	gl_draw::draw_arrow(
		um_bindings::um_vec(plane_center),
		um_bindings::um_vec(plane_center) + pq.normal() * 0.5,
		0.05, 8, 0.5, GEO::vec4f(1., 0.2, 0.1, 1.)
	);
}

void PathConstraintPaddingTool::hover_callback(double x, double y, int source) {
	if (/*step > 0 ||*/ !ctx.is_facet_hovered())
		return;

	facet_selector.paint(ctx.hovered_facet);
}

void PathConstraintPaddingTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void PathConstraintPaddingTool::scroll_callback(double xoffset, double yoffset) {}

void PathConstraintPaddingTool::validate_callback() {

	if (mode == Select) {

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