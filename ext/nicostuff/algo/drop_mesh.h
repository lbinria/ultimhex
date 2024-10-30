#ifndef DROP__MESH__H__
#define DROP__MESH__H__


#include <ultimaille/all.h>
#include <basic.h>

#include <algo/framework/trace.h>

#include <map>


using namespace UM;
#define ARG(Type,argument_name,defaut_value) \
	auto& _##argument_name(Type v) { __##argument_name = v; return *this; }; \
	Type __##argument_name = defaut_value;

#include <type_traits>
template<class M, class C>
struct Drop {
	Drop(M& m, C& attr) { std::cerr << "\nNot supposed to exists\n"; }
};

struct DropPointSet : public PointSetAttributes {
	DropPointSet(PointSet& pts, PointSetAttributes attr = {}) :pts(pts), PointSetAttributes(attr) {}

	PointSet& pts;

	ARG(std::string, just_save_filename, "");
	ARG(std::string, active_point_attribute, "");
	ARG(bool, lighting, true);


	template<class T>
	DropPointSet& add(PointAttribute<T>& c, std::string str) {
		points.push_back({ str,c.ptr });
		return *this;
	}

	DropPointSet& apply(std::string name = "") {
		if (!__just_save_filename.empty()) {
			write_by_extension(__just_save_filename, pts, *this);
			return *this;
		}
		if (!Trace::drop_mesh_is_active) return *this;
		std::string filename = Trace::outputdir + "/" +  name + "_" + std::to_string(Trace::num_drop++) + ".geogram";
		write_by_extension(filename, pts, *this);
		std::ofstream myfile;
		myfile.open(Trace::outputdir + "/view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 3'\n";
		myfile << "scene_graph.current().shader.colormap =  \"turbo;true;0;false;false;\"\n";

		if (!__active_point_attribute.empty()) {
			if (__active_point_attribute.compare("color") == 0) {
				myfile << "scene_graph.current().shader.painting = 'COLOR'\n";
				myfile << "scene_graph.current().shader.colors = 'vertices.color'\n";
				myfile << "scene_graph.current().shader.border_style = 'false'\n";
			}
			else {
				myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
				myfile << "scene_graph.current().shader.attribute = 'vertices." << __active_point_attribute << "'\n";
				myfile << "scene_graph.current().shader.autorange()\n";
			}
		}

		if (!__lighting)		myfile << "scene_graph.current().shader.lighting = 'false'\n";
		myfile.close();
		return *this;
	}
};



struct DropPolyLine : public PolyLineAttributes {
	DropPolyLine(PolyLine& pl, PolyLineAttributes attr = {}) :pl(pl), PolyLineAttributes(attr) {}
	PolyLine& pl;
	template<class T>
	DropPolyLine& add(PointAttribute<T>& c, std::string str) {
		points.push_back({ str,c.ptr });
		return *this;
	}
	template<class T>
	DropPolyLine& add(EdgeAttribute<T>& c, std::string str) {
		edges.push_back({ str,c.ptr });
		return *this;
	}


	ARG(std::string, active_point_attribute, "");
	ARG(std::string, active_segment_attribute, "");
	ARG(bool, lighting, true);
	ARG(std::string, just_save_filename, "");
	ARG(bool, show_vertices, true);

	DropPolyLine& apply(std::string name = "") {


		if (!__just_save_filename.empty()) {
			write_by_extension(__just_save_filename, pl, *this);
			return *this;
		}

		if (!Trace::drop_mesh_is_active) return *this;
		std::string filename = Trace::outputdir + "/" + name + "_" + std::to_string(Trace::num_drop++) + ".geogram";
		write_by_extension(filename, pl, *this);
		std::ofstream myfile;
		myfile.open(Trace::outputdir + "/view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		if (__show_vertices)
			myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 1'\n";
		myfile << "scene_graph.current().shader.colormap =  \"turbo;true;0;false;false;\"\n";

		if (!__active_point_attribute.empty()) {
			if (__active_point_attribute.compare("color") == 0) {
				myfile << "scene_graph.current().shader.painting = 'COLOR'\n";
				myfile << "scene_graph.current().shader.colors = 'vertices.color'\n";
			}
			else {
				myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
				myfile << "scene_graph.current().shader.attribute = 'vertices." << __active_point_attribute << "'\n";
				myfile << "scene_graph.current().shader.autorange()\n";
			}
		}
		if (!__active_segment_attribute.empty()) {
			if (__active_segment_attribute.compare("color") == 0) {
				myfile << "scene_graph.current().shader.painting = 'COLOR'\n";
				myfile << "scene_graph.current().shader.colors = 'edges.color'\n";
			}
			else {
				myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
				myfile << "scene_graph.current().shader.attribute = 'edges." << __active_segment_attribute << "'\n";
				myfile << "scene_graph.current().shader.autorange()\n";
			}
		}
		myfile.close();


		return *this;
	}
};



struct SurfaceAttributesBuilder : public SurfaceAttributes {
	SurfaceAttributesBuilder(SurfaceAttributes attr = {}) :SurfaceAttributes(attr) {}
	template<class T>
	SurfaceAttributesBuilder& add(FacetAttribute<T>& c, std::string str) {
		facets.push_back({ str,c.ptr });
		return *this;
	}
	template<class T>
	SurfaceAttributesBuilder& add(PointAttribute<T>& c, std::string str) {
		points.push_back({ str,c.ptr });
		return *this;
	}
	template<class T>
	SurfaceAttributesBuilder& add(CornerAttribute<T>& c, std::string str) {
		corners.push_back({ str,c.ptr });
		return *this;
	}
};

struct DropSurface : public SurfaceAttributes {
	DropSurface(Surface& m, SurfaceAttributes attr = {}) :m(m), SurfaceAttributes(attr) {}
	Surface& m;
	template<class T>
	DropSurface& add(FacetAttribute<T>& c, std::string str) {
		facets.push_back({ str,c.ptr });
		return *this;
	}
	template<class T>
	DropSurface& add(PointAttribute<T>& c, std::string str) {
		points.push_back({ str,c.ptr });
		return *this;
	}
	template<class T>
	DropSurface& add(CornerAttribute<T>& c, std::string str) {
		corners.push_back({ str,c.ptr });
		return *this;
	}




	ARG(std::string, active_facet_attribute, "");
	ARG(std::string, active_point_attribute, "");
	ARG(std::string, active_corner_attribute, "");
	ARG(std::string, active_uv_attribute, "");
	ARG(std::string, just_save_filename, "");

	ARG(bool, show_edges, true);
	ARG(bool, show_vertices, true);
	ARG(bool, show_border, true);
	ARG(bool, lighting, true);

	DropSurface& apply(std::string name = "") {


		if (!__just_save_filename.empty()) {
			write_by_extension(__just_save_filename, m, *this);
			return *this;
		}

		if (!Trace::drop_mesh_is_active) return *this;
		std::string filename = Trace::outputdir + "/" + name + "_" + std::to_string(Trace::num_drop++) + ".geogram";
		write_by_extension(filename, m, *this);
		std::ofstream myfile;
		myfile.open(Trace::lua_filename(), std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		myfile << "scene_graph.current().shader.colormap =  \"turbo;true;0;false;false;\"\n";

		if (!__active_uv_attribute.empty()) {
			myfile << "scene_graph.current().shader.painting = 'TEXTURE'\n";
			myfile << "scene_graph.current().shader.tex_coords = 'facet_corners." << __active_uv_attribute << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		else 		if (!__active_corner_attribute.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'facet_corners." << __active_corner_attribute << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		else if (!__active_facet_attribute.empty()) {
			if (__active_facet_attribute.compare("color") == 0) {
				myfile << "scene_graph.current().shader.painting = 'COLOR'\n";
				myfile << "scene_graph.current().shader.colors = 'facets.color'\n";
				myfile << "scene_graph.current().shader.border_style = 'false'\n";
			}
			else {
				myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
				myfile << "scene_graph.current().shader.attribute = 'facets." << __active_facet_attribute << "'\n";
				myfile << "scene_graph.current().shader.border_style = 'false'\n";
				myfile << "scene_graph.current().shader.autorange()\n";
			}
		}
		else if (!__active_point_attribute.empty()) {
			if (__active_point_attribute.compare("color") == 0) {
				myfile << "scene_graph.current().shader.painting = 'COLOR'\n";
				myfile << "scene_graph.current().shader.colors = 'vertices.color'\n";
				myfile << "scene_graph.current().shader.border_style = 'false'\n";
			}
			else {
				myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
				myfile << "scene_graph.current().shader.attribute = 'vertices." << __active_point_attribute << "'\n";
				myfile << "scene_graph.current().shader.autorange()\n";
			}
		}

		if (!__lighting)		myfile << "scene_graph.current().shader.lighting = 'false'\n";
		if (!__show_border)		myfile << "scene_graph.current().shader.border_style = 'true'\n";
		if (__show_vertices)	myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 3'\n";
		if (__show_edges)		myfile << "scene_graph.current().shader.mesh_style = 'true; .5 .5 .5 1; 1'\n";
		myfile.close();

		return *this;
	}
};




struct DropVolume : public VolumeAttributes {
	DropVolume(Volume& m) : m(m) {}
	DropVolume(Volume& m, VolumeAttributes attr) : VolumeAttributes(attr), m(m) {}

	template<class T> DropVolume& add(PointAttribute<T>& c, std::string str) {
		points.push_back({ str,c.ptr });
		return *this;
	}
	template<class T> DropVolume& add(CellAttribute<T>& c, std::string str) {
		cells.push_back({ str,c.ptr });
		return *this;
	}
	template<class T> DropVolume& add(CellFacetAttribute<T>& c, std::string str) {
		cell_facets.push_back({ str,c.ptr });
		return *this;
	}
	template<class T> DropVolume& add(CellCornerAttribute<T>& c, std::string str) {
		cell_corners.push_back({ str,c.ptr });
		return *this;
	}



	Volume& m;

	ARG(std::string, active_point_attribute, "");
	ARG(std::string, active_cell_attribute, "");
	ARG(std::string, active_cell_corner_attribute, "");
	ARG(std::string, just_save_filename, "");
	ARG(bool, show_vertices, false);

	DropVolume& apply(std::string name = "") {


		if (!__just_save_filename.empty()) {
			write_by_extension(__just_save_filename, m, *this);
			return *this;
		}

		if (!Trace::drop_mesh_is_active) return *this;
		std::string filename = Trace::outputdir + "/" + name + "_" + std::to_string(Trace::num_drop++) + ".geogram";
		write_by_extension(filename, m, *this);
		std::ofstream myfile;
		myfile.open(Trace::outputdir + "/view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 1'\n";
		myfile << "scene_graph.current().shader.colormap =  \"turbo;true;0;false;false;\"\n";
		if (!__active_cell_attribute.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'cells." << __active_cell_attribute << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		else if (!__active_point_attribute.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'vertices." << __active_point_attribute << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
			myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 3'\n";
		}
		else if (!__active_cell_corner_attribute.empty()) {
			myfile << "scene_graph.current().shader.painting = 'ATTRIBUTE'\n";
			myfile << "scene_graph.current().shader.attribute = 'cell_corners." << __active_cell_corner_attribute << "'\n";
			myfile << "scene_graph.current().shader.autorange()\n";
		}
		if (__show_vertices) myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 3'\n";
		else myfile << "scene_graph.current().shader.vertices_style = 'false; 0 0 0 1; 3'\n";
		myfile << "scene_graph.current().shader.shrink = 1\n";
		myfile.close();

		return *this;
	}
};






#endif