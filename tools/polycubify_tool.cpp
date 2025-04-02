#include "polycubify_tool.h"

#include <nicostuff/algo/framework/benjamin_API.h>

#include "../geom_ultimaille_binding.h"
#include "../context.h"
#include "../gl_draw.h"


void PolycubifyTool::run_nicostuff() {
	// Get UM cell facet attribute tet_flag from GEO mesh
	UM::CellFacetAttribute<int> tet_flag(ctx.tet, -1);
	um_bindings::um_attr_from_geo_attr<GEO::MESH_CELL_FACETS>(ctx.mesh_, "tet_flag", ctx.tet, tet_flag.ptr);

	// CellFacetAttribute<int> emb(ctx.hex, -1);

	try {
		// Polycubify and fill embedding
		ctx.emb_attr = std::make_unique<CellFacetAttribute<int>>(ctx.hex, -1);
		BenjaminAPI::polycubify(ctx.tet, tet_flag, ctx.hex, nhex_wanted, *ctx.emb_attr);

		// Make chart segmentation & init embedding
		ctx.tri_chart = std::make_unique<FacetAttribute<int>>(ctx.tet_bound->tri, -1);
		ctx.quad_chart = std::make_unique<FacetAttribute<int>>(ctx.hex_bound->quad, -1);
		BenjaminAPI::embeditinit(ctx.tet_bound->tri, *ctx.tri_chart, ctx.hex_bound->hex, *ctx.emb_attr, *ctx.quad_chart, false);

	} catch (const std::runtime_error &e) {
		Logger::warn("An error occur when trying to polycubify. Detail: " + std::string(e.what()));
		std::cout << "polycubify fail" << std::endl;
		return;
	}

	// Replace current GEO mesh by UM Hex
	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);

	// Replace current GEO mesh by UM Hex
	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);


	// Visualize feedback
	{
		FacetAttribute<int> surf_emb_attr(ctx.hex_bound->quad, -1);
		ctx.hex_bound->set_attribute_to_surface(*ctx.emb_attr, surf_emb_attr);
		write_by_extension("poly_emb.geogram", ctx.hex_bound->quad, {{}, {{"emb", surf_emb_attr.ptr}}, {}});
		write_by_extension("poly_trichart.geogram", ctx.tet_bound->tri, {{}, {{"trichart", ctx.tri_chart.get()->ptr}}, {}});
	}
	// TODO encapsulate in atomic unit ! + try catch to guarentee consistency

	// Write mesh

	// Save mesh metadata in json !!!!
	ctx.mesh_metadata = { 
		.filename = "polycubified.geogram", 
		.tet_filename = "polycubified.tet.geogram",
		.cell_type = GEO::MESH_HEX, 
		.attributes = {}
	};
	write_by_extension(ctx.mesh_metadata.filename, ctx.hex_bound->hex, {{}, {}, {{"emb", ctx.emb_attr->ptr}}, {}});
	write_by_extension(ctx.mesh_metadata.tet_filename, ctx.tet_bound->tet, {{}, {}, {{"chart", ctx.tri_chart->ptr}}, {}});
	ctx.mesh_metadata.save();

	// View
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	ctx.view.change_mode(ViewBinding::Mode::Volume);
}

void PolycubifyTool::run_robust_polycube() {

	CellFacetAttribute<int> flag(ctx.tet, -1);
	um_bindings::um_attr_from_geo_attr<MeshElementsFlags::MESH_CELL_FACETS>(ctx.mesh_, "tet_flag", ctx.tet, flag.ptr);

	// Write mesh file
	std::string filename = ctx.mesh_metadata.filename + ".mesh";
	
	write_by_extension(filename, ctx.tet, {{}, {}, {}, {}});

	std::string flag_filename = filename + ".flag";
	std::ofstream flag_file;
	flag_file.open(flag_filename);
	if (flag_file.is_open()) { 

		for (auto f : ctx.tet.iter_facets()) {
			// Convert flagging for robust polycube convention (see: https://github.com/fprotais/robustPolycube/tree/main)
			// int flag_map[] = {-1,1,0,3,2,5,4};
			int flag_map[] = {-1,1,3,5,0,2,4};
			int flag_val = flag_map[flag[f] + 1];

			flag_file << flag_val << "\n";
		}

		
		
		flag_file.close();
	} else {
		std::cout << "Failed to write flag file." << std::endl;
		return;
	}


	// Execute deformation
	std::string deform_cmd = std::string("/home/tex/Projects/robustPolycube/build/rb_generate_deformation") 
		+ " " + filename
		+ " " + flag_filename
		+ " " + filename + ".remesh.mesh"
		+ " " + filename + ".remesh.flag"
		+ " " + filename + ".polycuboid.mesh";


	// Execute quantization
	std::string quant_cmd = std::string("/home/tex/Projects/robustPolycube/build/rb_generate_quantization") 
		+ " " + filename + ".remesh.mesh"
		+ " " + filename + ".remesh.flag"
		+ " " + filename + ".polycuboid.mesh"
		+ " 1."
		+ " " + filename + ".hexmesh.mesh";


	// Execute post-processing
	std::string post_cmd = std::string("/home/tex/Projects/robustPolycube/build/rb_perform_postprocessing") 
		+ " " + filename + ".remesh.mesh"
		+ " " + filename + ".hexmesh.mesh"
		+ " " + filename + ".hexmesh_improved.mesh";

	int result = system(deform_cmd.c_str());

	if (result != 0) {
		std::cout << "Error on deformation" << std::endl;
		return;
	}

	result = system(quant_cmd.c_str());
	
	if (result != 0) {
		std::cout << "Error on quantization" << std::endl;
		return;
	}
	
	result = system(post_cmd.c_str());

	if (result != 0) {
		std::cout << "Error on post-processing" << std::endl;
		return;
	}

	// Load result
	
	// read_by_extension(filename + ".hexmesh_improved.mesh", ctx.hex);
	read_by_extension(filename + ".hexmesh.mesh", ctx.hex);
	ctx.hex.connect();

	// normalize
	BBox3 bb;
	for (auto v : ctx.hex.iter_vertices())
		bb.add(v.pos());
	
	UM::vec3 lbb = bb.max - bb.min;

	double max_coord = std::max(std::max(lbb.x, lbb.y), lbb.z);

	// normalize
	for (auto v : ctx.hex.iter_vertices()) {
		v.pos() -= bb.min;
		v.pos() /= max_coord;
	}

	// TODO below boilerplate code ! should refactor

	// Replace current GEO mesh by UM Hex
	ctx.hex_bound = std::make_unique<MyHexBoundary>(ctx.hex);

	um_bindings::geo_mesh_from_hexboundary(*ctx.hex_bound, ctx.mesh_);

	// TODO encapsulate in atomic unit ! + try catch to guarentee consistency

	// Write mesh

	// Save mesh metadata in json !!!!
	ctx.mesh_metadata = { 
		.filename = "polycubified.geogram", 
		.cell_type = GEO::MESH_HEX, 
		.attributes = {} 
	};

	write_by_extension(ctx.mesh_metadata.filename, ctx.hex, {{}, {}, {}, {}});
	ctx.mesh_metadata.save();

	// View
	ctx.mesh_gfx_.set_mesh(&ctx.mesh_);
	ctx.view.change_mode(ViewBinding::Mode::Volume);

}

bool PolycubifyTool::draw_object_properties() {

	// if (ImGui::CollapsingHeader("Polycubify")) {

		// Paint algo combo box selection
		ImGui::Text("Algorithm");
		if (ImGui::BeginCombo("Choose algo", algos[current_algo])) {

			for (int i = 0; i < IM_ARRAYSIZE(algos); i++)
			{
				bool isSelected = (current_algo == i);
				if (ImGui::Selectable(algos[i], isSelected))
				{
					current_algo = i;
				}

				if (isSelected)
					ImGui::SetItemDefaultFocus();
			}

			ImGui::EndCombo();
		}

		if (current_algo == Nicostuff) {

			ImGui::InputInt("Nb hex", &nhex_wanted);



		} else if (current_algo == RobustPolycube) {
			// Parameters here !
		}

		if (ImGui::BeginPopupModal("Mesh validity failed", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
			ImGui::TextUnformatted("This action requires all mesh faces to be flagged.");
			if (ImGui::Button("Close", ImVec2(120, 0))) {
				ImGui::CloseCurrentPopup();
			}
			ImGui::EndPopup();
		}

		if (ImGui::Button("Polycubify !")) {

			if (!check_flag_validity()) {
				ImGui::OpenPopup("Mesh validity failed");
				return false;
			}

			if (current_algo == PolycubeAlgo::Nicostuff)
				run_nicostuff();
			else if (current_algo == PolycubeAlgo::RobustPolycube) {
				run_robust_polycube();
			}

			return true;
		}

	// }

	return false;
}

bool PolycubifyTool::check_flag_validity() {
	// Check that all faces were flagged
	GEO::Attribute<GEO::signed_index_t> flag(
		ctx.mesh_.facets.attributes(), "flag"
	);

	for (auto f : ctx.mesh_.facets) {
		if (flag[f] < 0)
			return false;
	}

	return true;
}

void PolycubifyTool::draw_viewer_properties() {

}

void PolycubifyTool::draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) {

}

void PolycubifyTool::hover_callback(double x, double y, int source) {

}

void PolycubifyTool::mouse_button_callback(int button, int action, int mods, int source) {

}

void PolycubifyTool::scroll_callback(double xoffset, double yoffset) {}

void PolycubifyTool::validate_callback() {

}

void PolycubifyTool::escape_callback() {
	
}

bool PolycubifyTool::is_compatible() { 
	
	auto mesh_metadata_attr = ctx.mesh_metadata.get_attr("tet_flag");

	return ctx.mesh_metadata.cell_type == MESH_TET 
		&& mesh_metadata_attr.has_value() 
		&& mesh_metadata_attr.value().where == GEO::MESH_CELL_FACETS;
}