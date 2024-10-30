#ifndef HEX_STAT__H__
#define HEX_STAT__H__

#include <ultimaille/all.h>
#include <framework/trace.h>


#include "toolbox.h"
#include <iostream>

#include "mesh_geom.h"





struct HexRemeshQuality {
	HexRemeshQuality(Tetrahedra& tet, Hexahedra& hex, double eps_relative_to_tet_ave_edge_size = 1., bool show_max_dist = false);
	void drop_on_trace(std::string prefix = "algo");
	void drop_on_shell();
	double max_dist_to_tet;
	double max_dist_to_hex;
	double minSJ;
	double aveSJ;
	int edge_valence[10]; // 10+ goes to 9, boundary goes to 0
};


#endif