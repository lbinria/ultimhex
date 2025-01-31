#ifndef TOOLBOX_HEXAHEDRA_H__
#define TOOLBOX_HEXAHEDRA_H__


#include<algo/toolbox.h>
#include<algo/toolbox_pointset.h>

	struct ToolBoxHexahedra {
		ToolBoxHexahedra(Hexahedra& hex);
		double ave_edge_size();

		void copy_from(Hexahedra& other, bool duplicate_vertices = true);
		Quad3 facet_geom(int f);
		Hexahedron cell_geom(int c);
		void boundary(Quads& quad, bool shared_points);

		void eval_quality(CellAttribute<double>& SJ);
		void drop_scaled_jacobien(std::string name = "SJ");
		void drop_singular_edges(std::string name = "singu", bool show_hex = true);
		void drop_boundary(std::string name = "hex_boundary");



		void generate_regular_grid(int nb_subdivision);
		void generate_grid_with_singu_3(int nb_subdivision);
		void generate_grid_with_singu_5(int nb_subdivision);

		void split_all(int nb_iter = 1);

	/*
	* WARNING: does not support degenerated geometry
	*/
		void split();
		void read_best_efforts(std::string filename, VolumeAttributes& attribs);
		void read_best_efforts(std::string filename);
		void merge_vertices(DisjointSet& ds);

	Hexahedra& hex;
};


template<>
struct ToolBox<Hexahedra> : public ToolBoxHexahedra {
	ToolBox(Hexahedra& hex) :ToolBoxHexahedra(hex) {  }

};
#endif