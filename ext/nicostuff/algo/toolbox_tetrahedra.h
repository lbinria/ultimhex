#ifndef TOOL_BOX_TETRAHEDRA__H__
#define TOOL_BOX_TETRAHEDRA__H__
#include<algo/toolbox.h>
#include <drop_glyph.h>
#include <drop_attribute.h>
#include <toolbox_triangles.h>


struct ToolBoxTetrahedra {
	ToolBoxTetrahedra(Tetrahedra& tets);

	double ave_edge_size();

	void copy_from(Tetrahedra& other, bool connect = true);

	Triangle3 facet_geom(int f);

	double get_volume();


	vec3 facet_normal(int c, int lf);

	void boundary(Triangles& tri, bool shared_points);

	void read_best_efforts(std::string filename, VolumeAttributes& attribs, bool delete_isolated_vertices = false);
	void read_best_efforts(std::string filename, bool delete_isolated_vertices = false);

	void check_tet_is_manifold();

	void drop_tangled_cells(std::string name = "tangled");
	void drop_volume(std::string name = "tet vol");
	void drop_angles_quality(std::string name = "angle quality"	,double threshold= M_PI);

	void peel(double diheral_angle_treshold = M_PI / 16.);
	Tetrahedra& tets;
};
template<>
struct ToolBox<Tetrahedra> : public ToolBoxTetrahedra {
	ToolBox(Tetrahedra& tets) :ToolBoxTetrahedra(tets) {  }
};

#endif