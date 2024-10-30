#include "feature_curve_detector.h"

 FeatureEdgeDetector::FeatureEdgeDetector(Triangles& m) : m(m), hardness(m, 0) { drop = false; }

 FeatureEdgeDetector& FeatureEdgeDetector::verbose(bool b) {
	drop = b;
	return *this;
}

 FeatureEdgeDetector& FeatureEdgeDetector::dihedral_angle() {
	for (auto h : m.iter_halfedges()) {
		auto opp = h.opposite();
		if (!opp.active()) hardness[h] = 0;
		else {
			vec3 n0 = Triangle3(h.facet()).normal();
			vec3 n1 = Triangle3(opp.facet()).normal();
			hardness[h] = std::atan2(cross(n0, n1).norm(), n0 * n1) / M_PI;
		}
	}
	if (drop) Drop(m, hardness)._wireframe(true).apply_half_edge("dihedral angle");
	return *this;
}

 FeatureEdgeDetector& FeatureEdgeDetector::threshold(double angle) {
	for (auto h : m.iter_halfedges()) hardness[h] = (hardness[h] < angle) ? 0 : 1;
	if (drop) Drop(m, hardness)._skip_value(0)._wireframe(true).apply_half_edge("threshold");
	return *this;
}

 FeatureEdgeDetector& FeatureEdgeDetector::remove_small_features() {
	double ave = ToolBox(m).ave_edge_size();
	CornerAttribute<int> feature(m);
	int nfeats = apply(feature);
	std::vector<double> feat_length(nfeats);
	std::vector<int> ndandling_edge(nfeats);
	for (auto h : m.iter_halfedges()) if (feature[h] != -1) feat_length[feature[h]] += Segment3(h).vector().norm();
	for (auto h : m.iter_halfedges()) if (feature[h] != -1) {
		bool dandling = true;
		for (auto cir : h.to().iter_halfedges())
			dandling = dandling && (cir.to() == h.from() || feature[cir] == -1);
		if (dandling) ndandling_edge[feature[h]]++;
	}
	for (auto h : m.iter_halfedges()) if (feature[h] != -1) {
		if (ndandling_edge[feature[h]] == 0 && feat_length[feature[h]]<2. * ave) hardness[h] = 0;
		if (ndandling_edge[feature[h]]>0 && feat_length[feature[h]]<10. * ave) hardness[h] = 0;
	}
	if (drop) Drop(m, hardness)._skip_value(0)._wireframe(true).apply_half_edge("rm meaningless features");
	return *this;
}

 int FeatureEdgeDetector::apply(CornerAttribute<int>& feature, bool separate_halfedges) {
	DisjointSetWithNeutral ds(m.ncorners());
	for (auto h : m.iter_halfedges())
		if (hardness[h] < .5) ds.set_neutral(h);
		else if (!h.on_boundary() && !separate_halfedges) ds.merge(h, h.opposite());

	for (auto v : m.iter_vertices()) {
		int nfeats = 0;
		for (auto h0 : v.iter_halfedges()) if (hardness[h0] > .5) nfeats++;
		if (nfeats > 2) continue;

		for (auto h0 : v.iter_halfedges()) for (auto h1 : v.iter_halfedges()) {
			auto next = h1.prev();
			if (h0.opposite() == next) continue;
			if (hardness[h0] < .5) continue;
			if (hardness[next] < .5) continue;
			if (Geom3::vector_angle(Segment3(h0).vector(), Segment3(next).vector())>M_PI / 3.) continue;
			ds.merge(h0, next);
		}
	}

	return ds.get_sets_id(feature.ptr->data);
}

void test_FeatureEdgeDetector(Triangles& m) {
	CornerAttribute<int> feature(m);
	FeatureEdgeDetector(m).dihedral_angle().threshold().remove_small_features().remove_small_features().remove_small_features().apply(feature, false);
	//plop(nb_feats);
	Drop(m, feature)._wireframe(true).apply_half_edge("features");
}
