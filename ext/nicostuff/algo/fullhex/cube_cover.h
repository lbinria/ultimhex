#ifndef CUBE_COVER__H__
#define CUBE_COVER__H__
#include <ultimaille/all.h>
#include <fullhex/frame_field_3d.h>

struct CubeCover {
	CubeCover(Tetrahedra& p_m, FF3D& p_ff, CellCornerAttribute<vec3>& U,CellAttribute<bool>& locked_cells,CellFacetAttribute<int>& constraint_type);

	void cut_graph_3D_and_brush_ff(CellFacetAttribute<bool>& cut);
	void add_real_constraints(CellFacetAttribute<bool>& cut, ConstrainedLeastSquares& ls);
	void add_int_constraints(CellFacetAttribute<bool>& cut, ConstrainedLeastSquares& ls,CellCornerAttribute<vec3>& U_brush, bool force_boundary = false);
	void integrate_field( ConstrainedLeastSquares& ls, EdgeGraph& eg,double edge_length);

	void compute_GP();
	void compute_seamless();

	CellAttribute<bool>& locked_cells;
	Tetrahedra& m;
	CellFacetAttribute<int>& constraint_type;
	FF3D& ff;
	CellCornerAttribute<vec3>& U;
};

#endif