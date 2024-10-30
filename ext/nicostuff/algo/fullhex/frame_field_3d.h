#ifndef FRAME_FIELD_3D__H__
#define FRAME_FIELD_3D__H__
#include <ultimaille/all.h>
#include <volume/frame3D.h>


struct FF3D {
	FF3D(Tetrahedra& m, CellAttribute<mat3x3>& J) ;

	// basic access
	mat<3, 3>& operator[](int i) ;
	AxisPermutation permutation_around_edge(Volume::Halfedge h);
	bool edge_is_singular(Volume::Halfedge h) ;
	int singular_edge_stable_coordinate(Volume::Halfedge h) ;
	bool trace_streamline(int root_cell, int branch, vec3 P, std::vector<int>& cells, std::vector<vec3>& pts, int max_length = 100) ;
	void update_axis_permutation() ;

	// rendering
	void show_streamline(std::string name,int nb_shots = 1000) ;
	void show_singularity_graph(std::string name,bool with_neig_cubes = false) ;
	void show_cubes(std::string name,double scale = 1., bool with_flag= true) ;
	void show(std::string prefix= "") ;

	
	// computing 
	void init_with_constant_frame() ;
	void init_from_uvw(CellCornerAttribute<vec3>& U,bool smooth=true);


	
	void compute_with_fibers(CellAttribute<int>& fiber,CellFacetAttribute<int>& constraint_type,CellAttribute<bool> &locked,double data_fitting_w = 1.,int nb_pass = 1) ;
	bool apply_with_fibers(CellCornerAttribute<vec3>& U,CellAttribute<int>& fibers,CellFacetAttribute<int>& constraint_type,double data_fitting_w = 1.);	

	void compute_with_fibers_very_smooth(CellAttribute<int>& fiber, CellFacetAttribute<int>& constraint_type);
	bool apply_with_fibers_very_smooth(CellCornerAttribute<vec3>& U,CellAttribute<int>& fibers,CellFacetAttribute<int>& constraint_type);



	// post-process
	bool optimize_topology() ;
	void smooth_geom_only() ;
	


	
	void apply_ls(CellFacetAttribute<bool>& lock_normal,  double normal_coeff = 100.,CellAttribute<bool>* lock_cells = NULL);

	Tetrahedra& m;
	CellAttribute<mat3x3>& J;
	
	// ap defines the topology and is just a reference because it is not shared with other steps (it is reevaluated from J)
	CellFacetAttribute<AxisPermutation> ap;
};

#endif