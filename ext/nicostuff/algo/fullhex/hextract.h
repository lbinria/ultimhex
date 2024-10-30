#ifndef HEXTRACT__H__
#define HEXTRACT__H__
#include <ultimaille/all.h>
#include <create_shape.h>
#include <framework/trace.h>
#include "toolbox.h"
#include "drop_attribute.h"
#include "drop_glyph.h"
#include <fullhex/gp_basic.h>

#include <iostream>

#include <algo/volume/hex_edit.h>
#include <algo/mesh_geom.h>
#include<algo/volume/hex_stat.h>
#include<algo/dynamicKNN.h>

struct DualContour {
	DualContour(Hexahedra& hex) ;
	DualContour& init_hex_from_uvw(Tetrahedra& tet, CellCornerAttribute<vec3>& U) ;
	DualContour& fix_some_degenerated_cases(bool safe_mode=false);
	DualContour& compute_charts(HexBoundary& bound,FacetAttribute<int>& chart,int& nb_charts);
	bool try_to_pad();
	DualContour& pad_for_hex_edge() ;
	DualContour& pad_for_quad_corner() ;
	DualContour& postpad_propagate_constraints();
	
	DualContour& drop_SJ(std::string name ="");
	DualContour& show_constraints(std::string name);
	DualContour& untangle(double fit_scale = 1);
	// 0: no output, 1 only SJ, 2 all substeps
	DualContour& verbose_level(int i) { verbose = i; return *this; }

	std::tuple<int,vec3,mat3x3> get_fitting_constraint(std::vector<vec3>& pt,std::vector<vec3>& n);

	struct FacetConstraint {
		vec3 pt;
		vec3 n;
	};

	int verbose;
	Hexahedra& hex;
	CellFacetAttribute<FacetConstraint > cond;
	PointAttribute<bool> modified;
};



struct UntanglerHexTan {
	UntanglerHexTan(Hexahedra& m);
	void apply(DualContour& dc, double fit_scale = 1);
    Hexahedra& m;
};




#endif