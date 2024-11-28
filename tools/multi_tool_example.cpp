// #include "layer_pad_tool.h"

// #include <nicostuff/algo/framework/benjamin_API.h>

// #include "../geom_ultimaille_binding.h"
// #include "../context.h"
// #include "../gl_draw.h"

// #include "../helpers.h"


// std::string XTool::get_name() { return "X TOOL"; }

// bool XTool::draw_object_properties() {
// 	return false;
// }

// void XTool::draw_viewer_properties() {}

// void XTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

// }

// void XTool::hover_callback(double x, double y, int source) {

// }

// void XTool::mouse_button_callback(int button, int action, int mods, int source) {
// 	std::cout << "x tool ! " << std::endl;
// }

// void XTool::scroll_callback(double xoffset, double yoffset) {}

// void XTool::validate_callback() {

// }

// void XTool::key_callback(int key, int scancode, int action, int mods) {}

// void XTool::escape_callback() {

// }

// void XTool::clear() {

// }

// bool XTool::is_compatible() { 
// 	return true;
// }


// std::string LoopPadTool::get_name() { return "Layer padding"; }

// bool LoopPadTool::draw_object_properties() {
// 	return false;
// }

// void LoopPadTool::draw_viewer_properties() {}

// void LoopPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
// 	gl_draw::draw_path(hovered_path, hovered_color, true);
// 	gl_draw::draw_path(selected_path, selected_color, true);
// }

// void LoopPadTool::hover_callback(double x, double y, int source) {
// 	hovered_path.clear();

// 	if (ctx.hex.connected() && ctx.is_cell_hovered()) {

// 		Volume::Cell um_c(ctx.hex, ctx.hovered_cell);
		
// 		Volume::Halfedge hovered_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.hovered_edge, ctx.hovered_cell_lfacet)));
		
// 		if (hovered_he.active()) {

// 			helpers::loop_along(ctx.hex, hovered_he, [&](Volume::Halfedge &he, bool on_border) {
// 				UM::vec3 a = he.from().pos();
// 				UM::vec3 b = he.to().pos();
// 				hovered_path.push_back(a);
// 				hovered_path.push_back(b);
// 			});
// 		}
// 	}
// }

// void LoopPadTool::mouse_button_callback(int button, int action, int mods, int source) {
// 	selected_path = hovered_path;
// }

// void LoopPadTool::scroll_callback(double xoffset, double yoffset) {}

// void LoopPadTool::key_callback(int key, int scancode, int action, int mods) {}

// void LoopPadTool::validate_callback() {
// 	if (!ctx.is_cell_selected())
// 		return;

// 	CellFacetAttribute<bool> pad_face(ctx.hex);

// 	Volume::Cell um_c(ctx.hex, ctx.selected_cell);
// 	Volume::Halfedge start_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.selected_edge, ctx.selected_cell_lfacet)));

// 	// Measure time
// 	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

// 	helpers::layer_along(ctx.hex, start_he, [&](UM::Volume::Facet &f) {
// 		pad_face[f] = true;
// 	});

// 	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
// 	std::cout << "Loop cut duration = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;

// 	begin = std::chrono::steady_clock::now();

// 	BenjaminAPI::pad(ctx.hex, pad_face);

// 	end = std::chrono::steady_clock::now();

// 	std::cout << "Padding duration = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;


// 	um_bindings::geo_mesh_from_um_hex(ctx.hex, ctx.mesh_);
// 	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

// 	selected_path.clear();
// 	write_by_extension("padded.geogram", ctx.hex);
// }

// bool LoopPadTool::is_compatible() { 
// 	return ctx.mesh_metadata.cell_type == MESH_HEX;
// }

// void LoopPadTool::escape_callback() {
// 	clear();
// }

// void LoopPadTool::clear() {
// 	// Clear path
// 	hovered_path.clear();
// 	selected_path.clear();
// }





// bool LayerPadTool::draw_object_properties() {
// 	if (loop_pad_tool.is_compatible() && ImGui::Button("Loop padding")) {
// 		cur_tool = std::make_unique<LoopPadTool>(loop_pad_tool);
// 		cur_tool->clear();
// 		ctx.gui_mode = LayerPadding;
// 		return true;
// 	}

// 	if (xtool.is_compatible() && ImGui::Button("Zigzag padding")) {
// 		cur_tool = std::make_unique<XTool>(xtool);
// 		cur_tool->clear();
// 		ctx.gui_mode = LayerPadding;
// 		return true;
// 	}

// 	return false;
// }

// void LayerPadTool::draw_viewer_properties() {
// 	cur_tool->draw_viewer_properties();
// }

// void LayerPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
// 	if (cur_tool != NULL)
// 		cur_tool->draw(hovered_color, selected_color, colorMapInfo);
// }

// void LayerPadTool::hover_callback(double x, double y, int source) {
// 	if (cur_tool != NULL)
// 		cur_tool->hover_callback(x, y, source);
// }

// void LayerPadTool::mouse_button_callback(int button, int action, int mods, int source) {
// 	cur_tool->mouse_button_callback(button, action, mods, source);
// }

// void LayerPadTool::scroll_callback(double xoffset, double yoffset) {
// 	cur_tool->scroll_callback(xoffset, yoffset);
// }

// void LayerPadTool::validate_callback() {
// 	cur_tool->validate_callback();
// }

// bool LayerPadTool::is_compatible() { 
// 	return cur_tool->is_compatible();
// }

// void LayerPadTool::key_callback(int key, int scancode, int action, int mods) {
// 	cur_tool->key_callback(key, scancode, action, mods);
// }

// void LayerPadTool::escape_callback() {
// 	cur_tool->escape_callback();
// }

// void LayerPadTool::clear() {
// 	cur_tool->clear();
// }



// // bool LayerPadTool::draw_object_properties() {
// // 	if(ImGui::Button("Loop padding")) {
// // 		ctx.gui_mode = LayerPadding;
// // 		// TODO notify Change tool
// // 		hovered_path.clear();
// // 		selected_path.clear();
		
// // 		return true;

// // 	}

// // 	return false;
// // }

// // void LayerPadTool::draw_viewer_properties() {}

// // void LayerPadTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {
// // 	gl_draw::draw_path(hovered_path, hovered_color, true);
// // 	gl_draw::draw_path(selected_path, selected_color, true);
// // }

// // void LayerPadTool::hover_callback(double x, double y, int source) {
// // 	hovered_path.clear();

// // 	if (ctx.hex.connected() && ctx.is_cell_hovered()) {

// // 		Volume::Cell um_c(ctx.hex, ctx.hovered_cell);
		
// // 		Volume::Halfedge hovered_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.hovered_edge, ctx.hovered_cell_lfacet)));
		
// // 		if (hovered_he.active()) {

// // 			helpers::loop_along(ctx.hex, hovered_he, [&](Volume::Halfedge &he, bool on_border) {
// // 				UM::vec3 a = he.from().pos();
// // 				UM::vec3 b = he.to().pos();
// // 				hovered_path.push_back(a);
// // 				hovered_path.push_back(b);
// // 			});
// // 		}
// // 	}
// // }

// // void LayerPadTool::mouse_button_callback(int button, int action, int mods, int source) {
// // 	selected_path = hovered_path;
// // }

// // void LayerPadTool::scroll_callback(double xoffset, double yoffset) {}

// // void LayerPadTool::validate_callback() {

// // 	if (!ctx.is_cell_selected())
// // 		return;

// // 	CellFacetAttribute<bool> pad_face(ctx.hex);

// // 	Volume::Cell um_c(ctx.hex, ctx.selected_cell);
// // 	Volume::Halfedge start_he(ctx.hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(ctx.selected_edge, ctx.selected_cell_lfacet)));

// // 	// Measure time
// // 	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

// // 	helpers::layer_along(ctx.hex, start_he, [&](UM::Volume::Facet &f) {
// // 		pad_face[f] = true;
// // 	});

// // 	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
// // 	std::cout << "Loop cut duration = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;

// // 	begin = std::chrono::steady_clock::now();

// // 	BenjaminAPI::pad(ctx.hex, pad_face);

// // 	end = std::chrono::steady_clock::now();

// // 	std::cout << "Padding duration = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;


// // 	um_bindings::geo_mesh_from_um_hex(ctx.hex, ctx.mesh_);
// // 	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);

// // 	selected_path.clear();
// // 	write_by_extension("padded.geogram", ctx.hex);
// // }

// // bool LayerPadTool::is_compatible() { 
// // 	return ctx.mesh_metadata.cell_type == MESH_HEX;
// // }


// // Test extract layer
// // Quads q_out;
// // Volume::Cell um_c(hex, context_.hovered_cell);
// // Volume::Halfedge start_he(hex, um_c.halfedge(um_bindings::he_from_cell_e_lf(context_.selected_edge, context_.selected_lfacet)));
// // loop_cut(hex, start_he, [&](UM::Volume::Facet &f) {
	
// // 	int v_off = q_out.points.create_points(4);
// // 	q_out.points[v_off] = f.vertex(0).pos();
// // 	q_out.points[v_off + 1] = f.vertex(1).pos();
// // 	q_out.points[v_off + 2] = f.vertex(2).pos();
// // 	q_out.points[v_off + 3] = f.vertex(3).pos();

// // 	int f_off = q_out.create_facets(1);
// // 	q_out.vert(f_off, 0) = v_off;
// // 	q_out.vert(f_off, 1) = v_off + 1;
// // 	q_out.vert(f_off, 2) = v_off + 2;
// // 	q_out.vert(f_off, 3) = v_off + 3;

// // });
// // write_by_extension("result.geogram", q_out);