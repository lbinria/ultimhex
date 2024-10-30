#ifndef FEATURE_CURVE_DETECTOR__H_
#define FEATURE_CURVE_DETECTOR__H_

#include <fullhex/gp_basic.h> //OMG c'est quoi cette dependance ?!? ... c'est disjointSetWithNeutral, mais ca a rien a foutre la

struct FeatureEdgeDetector {
	FeatureEdgeDetector(Triangles& m);
	bool drop;
	FeatureEdgeDetector& verbose(bool b);

	FeatureEdgeDetector& dihedral_angle();

	FeatureEdgeDetector& threshold(double angle = .2);

	FeatureEdgeDetector& remove_small_features();

	int apply(CornerAttribute<int>& feature, bool separate_halfedges = false);
	CornerAttribute<double> hardness;
	Triangles& m;
};

void test_FeatureEdgeDetector(Triangles& m);

#endif