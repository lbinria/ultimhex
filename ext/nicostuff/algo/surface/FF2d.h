#ifndef __FF2D__H__
#define __FF2D__H__
#include <ultimaille/all.h>
#include <framework/trace.h>
#include "toolbox.h"
#include "drop_attribute.h"


struct FF2D {
	FF2D(Triangles& m, CornerAttribute<bool>& feature, CornerAttribute<vec2>& U, CornerAttribute<int>& LS) : m(m), feature(feature), U(U), LS(LS) {}

	// input/output
	Triangles& m;

	// output
	CornerAttribute<bool>& feature;
	CornerAttribute<vec2>& U;
	CornerAttribute<int>& LS;

	void show_grad(bool with_inverse = true);
	void show_iso(bool with_inverse = true);
	void show_vec2(CornerAttribute<vec2>& uv, FacetAttribute<vec2>& in_vec2, double scale = .1);
	void show_features();
	void show_singularities();
	void show();

	void detect_feature_by_dihedral_angle_threshold(double threshold = M_PI / 5.);
	void estimate_curvature(int nb_smoothing_iters = 5);
	void init_iso_map(CornerAttribute<vec2>& uv);
	void compute_LS(CornerAttribute<vec2>& metric);

	int index(Surface::Vertex v);

	void extrapolate_FF_from_feature(CornerAttribute<vec2>& metric);
	void smooth_FF(CornerAttribute<vec2>& metric);

	void init_map_for_sharp_corners(CornerAttribute<vec2>& iso_metric, CornerAttribute<vec2>& corner_metric);

	void use_case_CAD();
	void use_case_CAD_feature_is_given();
	void use_case_scanned_mesh();
};
#endif