#ifndef SEAMLESS_TO_GP__H__
#define SEAMLESS_TO_GP__H__
#include <ultimaille/all.h>

//
using namespace UM::Linear;

#include <fullhex/frame_field_3d.h>

struct SeamlessToGP {
	SeamlessToGP(Tetrahedra& p_m, CellCornerAttribute<vec3>& U,CellFacetAttribute<int>& constraint_type);

	void cut_graph_3D_and_brush_ff();
	void add_real_constraints(ConstrainedLeastSquares& ls, bool force_boundary = false);
	void add_int_constraints(ConstrainedLeastSquares& ls,bool force_boundary = false);
	void integrate_field( ConstrainedLeastSquares& ls, EdgeGraph& eg);

	void compute(bool force_boundary=true);
	
	Tetrahedra& m;
	CellFacetAttribute<bool> cut;
	CellFacetAttribute<int>& constraint_type;
	CellCornerAttribute<vec3>& U;
};

#endif