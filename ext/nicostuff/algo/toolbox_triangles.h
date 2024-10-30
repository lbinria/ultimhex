#ifndef TOOLBOX__TRIANGLE__H__
#define TOOLBOX__TRIANGLE__H__

#include <framework/trace.h>




	
struct ToolBoxTriangles {
	ToolBoxTriangles(Triangles& m) :m(m) {  }
	Triangles& m;

	void copy_from(Triangles& other, bool duplicate_vertices=true);

	double ave_edge_size();
	vec3 barycenter(int f);
	void split_triangle_4();
	void tetgen(Tetrahedra& tet, bool with_inner_vertices = false);
	void read_best_efforts(std::string filename, SurfaceAttributes& attribs);
	void read_best_efforts(std::string filename);
	bool make_orientable_manifold();
	int add_triangle(Triangle3 tr);
	void add_triangles(Triangles& m1);
	void merge_vertices(double epsilon = 0);
	void kill_degenerated_facets_topo();
	void kill_opposite_facets();

	struct TriangleConnectivityStat {
		TriangleConnectivityStat(Surface& m);
		std::string get_warning_string() ;
		std::string get_string() ;
		int nb_vertices;
		int nb_isolated_vertices;
		int nb_facets;
		int nb_samosas;
		int nb_null_edges;
		std::vector<int> multiple_ombrella_vertices;
		std::map<int, int > nb_opposites;
		std::map<int, int > nb_duplicated_edges;
	};

	TriangleConnectivityStat get_connectivity_info() {
		return TriangleConnectivityStat(m);
	}
	void laplasmooth(double coeff = .95, bool lock_non_manifold_edges = true);
};


template<>
struct ToolBox<Triangles> : public ToolBoxTriangles {
	ToolBox(Triangles& m) :ToolBoxTriangles(m) {  }
};

#endif