#ifndef TOOLBOX_POLYLINE__H__
#define TOOLBOX_POLYLINE__H__
#include<algo/toolbox.h>



	struct ToolBoxPolyLine {
		ToolBoxPolyLine(PolyLine& pl);
		void copy_from(PolyLine& other, bool duplicate_vertices = true) {
			if (duplicate_vertices) ToolBoxPointSet(pl.points).copy_from(other.points);
			else					pl.points = other.points;
			pl.create_edges(other.nedges());
			FOR(e, pl.nedges()) FOR(lv, 2) pl.vert(e, lv) = other.vert(e, lv);
		}
		void triangulate_with_shewchuck(Triangles& m, int verbose = false, int may_add_vertices = false);
		double ave_edge_size();
		void kill_isolated_vertices();
		int add_segment(vec3 A, vec3 B);
	PolyLine& pl;
};

	template<>
	struct ToolBox<PolyLine> : public ToolBoxPolyLine {
		ToolBox(PolyLine& pl) :ToolBoxPolyLine(pl) {  }
	};
#endif