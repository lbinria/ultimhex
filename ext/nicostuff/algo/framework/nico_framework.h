#ifndef NICO_FRAMEWORK__H__
#define NICO_FRAMEWORK__H__
#include <framework/param_parser.h>
#include <framework/meta.h>
#include "drop_attribute.h"
#include <regex>





struct NicoFramework;
void framework_main(NicoFramework &param);
void framework_parameters(NicoFramework &fw);
struct NicoFramework :Parameters{

    NicoFramework(int argc, char** argv,std::string input_file=""){
	    add("file", "output_name", "out.geogram").type_of_param("system").description("output name");
	    add("directory", "output_path", "C:/NICO/prog/output/tmp").type_of_param("system").description("the directory to ouputs: .lua, .csv and results");
	    add("directory", "result_path", "C:/NICO/prog/output/result").type_of_param("system").description("the directory to ouput results");
	    add("bool", "show_dropped_mesh", "true").type_of_param("system").description("Lauch a viewer to see dropped meshes");
		add("file", "graphite", "C:/NICO/prog/GraphiteThree/build/Windows/bin/Release/graphite.exe").type_of_param("system").description("Viewer graphite.exe location");
		add("file", "gmsh", "C:/NICO/prog/gmsh/gmsh.exe").type_of_param("system").description("gmsh.exe location");
		add("input", "input",input_file).type_of_param("system").description("The input mesh");
		framework_parameters(*this);
        init_from_args(argc, argv);
		Trace::initialize((*this)["output_path"], "");
		Trace::SwitchDropInScope::force_activity((*this)["show_dropped_mesh"] && !has_run_from());
    }

    Param&  add_parameter(std::string type, std::string name, std::string default_value){
        return add(type, name, default_value);
    }
	void set_default_input(std::string input_file){
		if (std::string(data["input"]).empty()) 
			data["input"].set(input_file);
	}

    ~NicoFramework(){
    	Trace::show_log();
	    Trace::append_log(std::string((*this)["output_path"]) + "/log.txt");
	    if (!std::filesystem::exists(Trace::lua_filename())) return;
		auto cmd = std::string((*this)["graphite"])+" " + Trace::lua_filename();
	    if ((*this)["show_dropped_mesh"]) system(cmd.c_str());
    }

	std::string output_file(){return std::string((*this)["result_path"])+"/"+std::string((*this)["output_name"]);}
};



void load_triangles(Triangles& m, SurfaceAttributes& attribs, NicoFramework& fw) {
	std::string extension = fw["input"];
	extension = extension.substr(extension.find_last_of(".") + 1);
	if (extension.compare("step") == 0) {
		std::ofstream geofile("tmp.geo"); // DIRTY: it will be created in working directory :(
		if (geofile.is_open()) {
			geofile << "Merge \"" << std::string(fw["input"]) << "\";\nMesh.Algorithm3D = 1;\nMesh  2;\nMesh.Format = 30;\nMesh.SaveAll = 1;\nSave \"tmp.mesh\";\n";
			geofile.close();
		}
		else Trace::abort("Failed to write the .geo script to convert .step into .mesh");
		std::string cmd = std::string(fw["gmsh"]) + " tmp.geo -clscale .15 -parse_and_exit";
		system(cmd.c_str());
		ToolBox(m).read_best_efforts("tmp.mesh", attribs);
	}
	else ToolBox(m).read_best_efforts(fw["input"], attribs);
}

void load_tetrahedra(Tetrahedra& m, VolumeAttributes& attribs, NicoFramework& fw) {
	std::string extension = fw["input"];
	extension = extension.substr(extension.find_last_of(".") + 1);
	if (extension.compare("step") == 0) {
		std::ofstream geofile("tmp.geo"); // DIRTY: it will be created in working directory :(
		if (geofile.is_open()) {
			geofile << "Merge \"" << std::string(fw["input"]) << "\";\nMesh.Algorithm3D = 1;\nMesh  3;\nMesh.Format = 30;\nMesh.SaveAll = 1;\nSave \"tmp.mesh\";\n";
			geofile.close();
		}
		else Trace::abort("Failed to write the .geo script to convert .step into .mesh");
		std::string cmd = std::string(fw["gmsh"]) + " tmp.geo -clscale .15 -parse_and_exit";
		system(cmd.c_str());
		ToolBox(m).read_best_efforts("tmp.mesh", attribs);
	}
	else ToolBox(m).read_best_efforts(fw["input"], attribs);
}


int main(int argc, char** argv) {
	std::cerr<<"--------------------------------PROGRAM-----------------------------------------------\n"<< argv[0]<<"\n";
	NicoFramework fw(argc, argv);
	try {
		std::cerr << "--------------------------------ARGUMENTS----------------------------------------------\n";
		std::cerr << std::regex_replace(fw.str_values(), std::regex(" "), "\n");
		std::cerr << "---------------------------------------------------------------------------------------\n\n";
		framework_main(fw);
	}
	catch (const std::exception& e) {
		Trace::log_string("END", e.what(), 3);
	}
	std::cerr<<"FINISH TEST_EXE\n";
}
#endif