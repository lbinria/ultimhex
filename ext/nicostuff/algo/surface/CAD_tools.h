#ifndef __CAD_TOOLS__H__
#define __CAD_TOOLS__H__
#include <ultimaille/all.h>


#include <surface/pointset_in_surface.h>


struct EdgeFlip {
	EdgeFlip(Triangles& m, Surface::Halfedge h);
	bool valid();
	void apply(CornerAttribute<bool>& feature);

	Triangles& m;
	Surface::Halfedge h;
	Surface::Halfedge opp;
};



struct EdgeCollapse {
	EdgeCollapse(Triangles& m, Surface::Halfedge h);
	bool valid();
	void apply(CornerAttribute<bool>& feature);
	void apply(FacetAttribute<int>& chart);
	Triangles& m;
	Surface::Halfedge h;
	Surface::Halfedge opp;
};




struct TrianglesOptim {
	TrianglesOptim(Triangles& tri, CornerAttribute<bool>& feature);

	void init_embedding();
	void smooth_pass();
	void edgeflip_pass();
	void edgecollapse_pass(bool feature_preserving = true);
	void apply();

	Triangles& tri;
	PointSetEmbedding emb;
	CornerAttribute<bool>& is_feature;
};


#endif