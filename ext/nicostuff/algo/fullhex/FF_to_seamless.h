#ifndef FF_TO_SEAMLESS__H__
#define FF_TO_SEAMLESS__H__
#include <ultimaille/all.h>


#include <fullhex/frame_field_3d.h>

struct FFToSeamless {
	FFToSeamless(Tetrahedra& p_m, FF3D& p_ff, CellCornerAttribute<vec3>& U,CellFacetAttribute<int>& constraint_type);

	void cut_graph_3D_and_brush_ff(CellFacetAttribute<bool>& cut);
	bool add_real_constraints(CellFacetAttribute<bool>& cut, ConstrainedLeastSquares& ls, bool force_boundary = false);
	bool integrate_field( ConstrainedLeastSquares& ls, double edge_length);

	bool compute(bool force_boundary=true);

	Tetrahedra& m;
	CellFacetAttribute<int>& constraint_type;
	FF3D& ff;
	CellCornerAttribute<vec3>& U;
};

#endif