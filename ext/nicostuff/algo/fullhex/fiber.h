#ifndef fiber__H__
#define fiber__H__

#include <ultimaille/all.h>
#include <volume/frame3D.h>
//#include <geom3.h>
#include <fullhex/gp_basic.h>

struct FiberSegmentation {
	FiberSegmentation(Tetrahedra& m, CellCornerAttribute<vec3>& U, CellAttribute<int>& fiber,CellAttribute<bool>& selection,CellFacetAttribute<int>& constraint_type) : m(m), U(U), fiber(fiber),selection(selection),bound(m),constraint_type(constraint_type) {}
	
	void compute_chart_id(FacetAttribute<int>& chart_id,FacetAttribute<vec3>& stable_direction,int &nb_charts) ;
	
	void compute_singu_constraints(CellAttribute<std::array<bool,3> > &constraint);
	bool grow_from_seed(CellAttribute<std::array<bool,3> >& constraint,CellAttribute<std::array<bool,3> >& in_previous, Volume::Facet seed_f_tet,CellAttribute<std::array<bool,3> >& already_in,CellFacetAttribute<int>& local_constraint_type, FacetAttribute<int>& chart_id,FacetAttribute<vec3>& stable_direction,int &nb_charts);


	void cleanup_fibers();

	void grow_regions(CellAttribute<std::array<bool,3> >& constraint);
	bool orient_seed(CellAttribute<std::array<bool,3> >& constraint,CellAttribute<std::array<bool,3> >& in_previous,Volume::Facet seed_f_tet);
	void trace_fiber();
	bool apply() ;
	bool apply_to_deformed_only();

	Tetrahedra& m;
	CellCornerAttribute<vec3>& U;
	CellAttribute<int>& fiber;
	CellAttribute<bool>& selection;
		/*
		* constraint_type -2 => locked 
		* constraint_type -1 => free 
		* constraint_type i \in 0,1,2 => grad(U[i]) is the normal 
		* constraint_type i \in 3,4,5 => grad(U[i-3]) is the stable direction
		*/
	CellFacetAttribute<int>& constraint_type;
	TetBoundary bound;

};






#endif