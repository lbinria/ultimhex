
#include <surface/CAD_tools.h>
#include <framework/trace.h>
#include "toolbox.h"
#include "drop_attribute.h"
 TrianglesOptim::TrianglesOptim(Triangles& tri, CornerAttribute<bool>& is_feature) : tri(tri), emb(tri.points), is_feature(is_feature) {
}

inline void TrianglesOptim::init_embedding() {
	{// init points embedding
		emb.init_from_triangles(tri, &is_feature);

		PointAttribute<int> nfeats(tri, 0);

		for (auto e : emb.pl_emb.iter_edges()) nfeats[e.from()]++;

		for (auto h : tri.iter_halfedges()) {
			Surface::Vertex v = h.from();

			// v is already by another halfedge
			if (emb.dim[v] != -1)
				continue;

			// v is dandling or intersection
			if (nfeats[v] == 1 || nfeats[v] > 2) {
				emb.dim[v] = 0; emb.id[v] = v;
			}
			// v doesn't touch any is_feature
			if (nfeats[v] == 0) {
				emb.dim[v] = 2; emb.id[v] = h.facet();
			}
		}

		// v is on a is_feature
		for (auto e : emb.pl_emb.iter_edges()) {
			auto v = e.from();
			int val = 0;
			for (auto cir : v.iter_edges()) val++;
			if (val == 2) { emb.dim[v] = 1; emb.id[v] = e; }
		}

		Drop(tri, nfeats).apply("nfeats");
	}
}

inline void TrianglesOptim::smooth_pass() {
	PointAttribute<vec3> smoothed_pos(emb.pts);
	{
		LeastSquares ls(3 * tri.nverts());
		for (auto v : tri.iter_vertices()) {
			auto constraints = emb.constrained_direction(v);
			for (vec3 c : constraints) {
				if (c.norm2() < 1e-20) continue;
				ls.add_to_energy(100. * (c[0] * X(3 * v + 0) + c[1] * X(3 * v + 1) + c[2] * X(3 * v + 2) - c * v.pos()));
			}
		}
		for (auto h : tri.iter_halfedges()) FOR(d, 3)
			ls.add_to_energy(X(3 * h.from() + d) - X(3 * h.to() + d));
		ls.solve();
		for (auto v : tri.iter_vertices()) FOR(d, 3) smoothed_pos[v][d] = ls.X[3 * v + d];
	}
	emb.move_toward(smoothed_pos);
	DropSurface(tri).apply("advect");
}

inline void TrianglesOptim::edgeflip_pass() {

	auto corner_angle = [&](Triangle3 t, int corner_id) {
		return Geom3::vector_angle(t[(corner_id + 1) % 3] - t[corner_id], t[(corner_id + 2) % 3] - t[corner_id]);
		};

	bool all_done = false;
	while (!all_done) {
		all_done = true;
		FOR(h_id, tri.ncorners()) {
			Surface::Halfedge h(tri, h_id);
			if (!h.active()) continue;
			if (is_feature[h]) continue;

			double flip_quality = 0;
			if (1) {// improve triangle corners min
				vec3 P[4] = { h.from().pos(),h.opposite().prev().from().pos(),h.to().pos(),h.prev().from().pos() };
				Triangle3 tr_avt[2] = { Triangle3(P[0],P[1],P[2]), Triangle3(P[0],P[2],P[3]) };
				Triangle3 tr_ap[2] = { Triangle3(P[0],P[1],P[3]), Triangle3(P[1],P[2],P[3]) };
				double min_angle_avt = M_PI;
				FOR(t, 2)FOR(lv, 3) min_angle_avt = std::min(min_angle_avt, corner_angle(tr_avt[t], lv));
				double min_angle_ap = M_PI;
				FOR(t, 2)FOR(lv, 3) min_angle_ap = std::min(min_angle_ap, corner_angle(tr_ap[t], lv));
				if (min_angle_ap > min_angle_avt
					&& tr_avt[0].normal() * tr_avt[1].normal() > 0
					&& tr_ap[0].normal() * tr_ap[1].normal() > 0
					)flip_quality = 1;
			}

			if (flip_quality <= 0) continue;
			EdgeFlip ef(tri, h);
			if (ef.valid()) {
				ef.apply(is_feature);
				all_done = false;
			}
		}
	}
	tri.compact();
	//Drop(tri, is_feature)._skip_value(false)._wireframe(true).apply_half_edge("is_feature");
	DropSurface(tri).apply("edgeflip");
}

inline void TrianglesOptim::edgecollapse_pass(bool is_feature_preserving ) {

	double ave = ToolBox(tri).ave_edge_size();
	if (is_feature_preserving) {
		for (double min_edge_size = ave / 20.; min_edge_size < ave / 2.; min_edge_size *= 2.) {
			bool all_done = false;

			while (!all_done) {
				all_done = true;
				FOR(h_id, tri.ncorners()) {
					Surface::Halfedge h(tri, h_id);
					if (!h.active()) continue;

					// do not produce incoherent is_feature
					if (is_feature[h.next()] != is_feature[h.prev()]) continue;
					if (is_feature[h.opposite().next()] != is_feature[h.opposite().prev()]) continue;

					// do not move is_feature
					if (!is_feature[h]) {
						bool we_have_a_problem = false;;
						for (auto cir : h.from().iter_halfedges()) if (cir != h && is_feature[cir])we_have_a_problem = true;
						if (we_have_a_problem) continue;
					}

					if (min_edge_size < (h.to().pos() - h.from().pos()).norm()) continue;

					EdgeCollapse ec(tri, h);
					if (ec.valid()) {
						ec.apply(is_feature);
						all_done = false;
					}
				}
			}
		}
		tri.compact();

	}
	else {

		FacetAttribute<int> chart(tri);
		DisjointSet ds(tri.nfacets());
		for (auto h : tri.iter_halfedges()) if (!is_feature[h]) ds.merge(h.facet(), h.opposite().facet());
		ds.get_sets_id(chart.ptr->data);
		//Drop(tri, chart).apply("charts");



		auto collapse_cost = [&](Surface::Halfedge h) {
			EdgeCollapse ec(tri, h);
			if (!ec.valid()) return 1e20;
			double sum = 0;
			vec3 np = .5 * (h.to().pos() + h.from().pos());
			for (auto cir : h.from().iter_halfedges()) {
				if (h == cir) continue;
				if (h.opposite() == cir.prev()) continue;
				Tetrahedron tet(np, h.to().pos(), cir.to().pos(), cir.next().to().pos());
				sum += std::abs(tet.volume());
			}
			for (auto cir : h.to().iter_halfedges()) {
				if (h == cir.prev()) continue;
				if (h.opposite() == cir) continue;
				Tetrahedron tet(np, h.from().pos(), cir.to().pos(), cir.next().to().pos());
				sum += std::abs(tet.volume());
			}

			return sum * Segment3(h).length();
			};
		//auto collapse_cost = [&](Surface::Halfedge h) {
		//	EdgeCollapse ec(tri, h);
		//	if (!ec.valid()) return 1e20;
		//	double sum = 0;
		//	for (auto cir : h.from().iter_halfedges()) {
		//		if (h == cir) continue;
		//		if (h.opposite() == cir.prev()) continue;
		//		Tetrahedron tet(h.from().pos(),h.to().pos(),cir.to().pos(), cir.next().to().pos());
		//		sum += std::abs(tet.volume());
		//	}

		//	return sum*Segment3(h).length();
		//};
		CornerAttribute<double> stored_collapse_cost(tri, 1e20);
		for (auto h : tri.iter_halfedges()) stored_collapse_cost[h] = collapse_cost(h);
		//Drop(tri, stored_collapse_cost).apply_interpolated("eccost");

		while (true) {
			//static int it = 0; if (it++ % 100 == 0) {tri.compact();Drop(tri, chart).apply("charts");}

			Surface::Halfedge h(tri, -1); 
			{// find min _cost h
				double min_cost = 1e20;
				for (auto hh : tri.iter_halfedges()) if (min_cost > stored_collapse_cost[hh]) {
					min_cost = stored_collapse_cost[hh];
					h = hh;
				}
				if (min_cost > 400) break;							
			}
			EdgeCollapse ec(tri, h);
			um_assert(ec.valid());
			int fid = tri.nfacets();
			// opt middle
			h.to().pos() = .5 * (h.to().pos() + h.from().pos());
			auto v = h.to();
			ec.apply(chart);

			for (auto h : v.iter_halfedges()) {
				for (auto cir : h.to().iter_halfedges()) {
					stored_collapse_cost[cir] = collapse_cost(cir);
					stored_collapse_cost[cir.opposite()] = collapse_cost(cir.opposite());
				}
			}

			//Surface::Facet f(tri, fid);
			//while (f<tri.nfacets()) {
			//	for (auto h : f.iter_halfedges()) {
			//		for (auto cir : h.from().iter_halfedges())
			//			stored_collapse_cost[cir] = collapse_cost(cir);
			//	}
			// 
			//	f++;
			//}
		}
		tri.compact();

		CornerAttribute<int> feature(tri);
		FeatureEdgeDetector(tri).dihedral_angle().threshold(.1)
			//.remove_small_features().remove_small_features().remove_small_features()
			.apply(feature, false);
		Drop(tri,feature).apply("featurecurves");
		for (auto h : tri.iter_halfedges()) is_feature[h] = 
			//(feature[h] != -1);
			(chart[h.facet()] != chart[h.opposite().facet()]);
		//Drop(tri, chart).apply("charts");
	}
}

 void TrianglesOptim::apply() {
	init_embedding();

	DropSurface(tri).apply("trin");
	FOR(it, 4) {
		plop(it);
		smooth_pass();
		edgeflip_pass();
		edgecollapse_pass();
	}
	DropSurface(tri).apply("tri");

}

inline EdgeFlip::EdgeFlip(Triangles& m, Surface::Halfedge h) : m(m), h(h), opp(h.opposite()) { }

inline bool EdgeFlip::valid() {
	Surface::Vertex from = h.from();
	Surface::Vertex to = h.to();
	if (!opp.active()) return false;
	if (ToolBox(to).valence() < 4) return false;
	if (ToolBox(from).valence() < 4) return false;
	for (auto cir_left : h.next().to().iter_halfedges())
		if (cir_left.to() == opp.next().to())
			return false;;
	return true;
}

inline void EdgeFlip::apply(CornerAttribute<bool>& is_feature) {
	int old_f[2] = { h.facet(),opp.facet() };

	// create facets
	auto nf = m.Surface::conn->create_facet({ h.to(),h.prev().from(),opp.prev().from() });
	auto nopp = m.Surface::conn->create_facet({ opp.to(),opp.prev().from(),h.prev().from() });

	// transfert attribute
	is_feature[nf.halfedge(0)] = is_feature[h.next()];
	is_feature[nf.halfedge(2)] = is_feature[opp.prev()];

	is_feature[nopp.halfedge(0)] = is_feature[opp.next()];
	is_feature[nopp.halfedge(2)] = is_feature[h.prev()];

	// deactivate old facets
	FOR(f, 2) m.Surface::conn->active[old_f[f]] = false;
}

inline EdgeCollapse::EdgeCollapse(Triangles& m, Surface::Halfedge h) :m(m), h(h), opp(h.opposite()) { }

inline bool EdgeCollapse::valid() {
	if (opp.active())	for (auto cir_from : h.from().iter_halfedges()) if (!cir_from.opposite().active()) return false;
	for (auto cir_from : h.from().iter_halfedges())
		for (auto cir_to : h.to().iter_halfedges())
			if (cir_from.to() == cir_to.to())
				if (cir_from.to() != h.next().to())
					if (!opp.active() || cir_from.to() != opp.next().to())
						return false;
	return true;
}

inline void EdgeCollapse::apply(CornerAttribute<bool>& is_feature) {
	std::vector<int> h_incident;
	for (auto cir_from : h.from().iter_halfedges())
		h_incident.push_back(cir_from);


	for (int it_id : h_incident) {
		Surface::Halfedge it(m, it_id);

		if (it.to() == h.to()) continue;
		if (it.prev().from() == h.to()) continue;
		auto nf = m.Surface::conn->create_facet({ h.to(), it.to(), it.prev().from() });
		is_feature[nf.halfedge(0)] = is_feature[it];
		is_feature[nf.halfedge(1)] = is_feature[it.next()];
		is_feature[nf.halfedge(2)] = is_feature[it.prev()];
	}
	for (int it_id : h_incident)
		m.Surface::conn->active[Surface::Halfedge(m, it_id).facet()] = false;
}
inline void EdgeCollapse::apply(FacetAttribute<int>& chart) {
	std::vector<int> h_incident;
	for (auto cir_from : h.from().iter_halfedges())
		h_incident.push_back(cir_from);


	for (int it_id : h_incident) {
		Surface::Halfedge it(m, it_id);
		if (it.to() == h.to()) continue;
		if (it.prev().from() == h.to()) continue;
		auto nf = m.Surface::conn->create_facet({ h.to(), it.to(), it.prev().from() });
		chart[nf] = chart[it.facet()];
	}
	for (int it_id : h_incident)
		m.Surface::conn->active[Surface::Halfedge(m, it_id).facet()] = false;
}


