#include<algo/volume/hex_stat.h>
#include<algo/drop_attribute.h>
#include<algo/drop_glyph.h>


HexRemeshQuality::HexRemeshQuality(Tetrahedra& tet, Hexahedra& hex, double eps_relative_to_tet_ave_edge_size, bool show_max_dist  ) {
	if (tet.ncells()==0 || hex.ncells()==0) return;
	hex.connect();
	//return;
	{// eval min Scaled Jacobian
		EdgeGraph eg(hex);
	
		EdgeAttribute<int> valence(eg, 0);

		for (auto h : hex.iter_halfedges())  {
			auto e = eg.edge_from_halfedge(h);
			valence[eg.edge_from_halfedge(h)]++;
		}
		
		for (auto h : hex.iter_halfedges())  if (!h.facet().opposite().active()) valence[eg.edge_from_halfedge(h)] = 0;

		FOR(i, 10) edge_valence[i] = 0;
		for (auto e : eg.iter_edges()) edge_valence[std::min(9, valence[e])]++;
	}


		double eps = ToolBox(tet).ave_edge_size() * eps_relative_to_tet_ave_edge_size;

		{// eval min Scaled Jacobian
			CellAttribute<double> SJ(hex);
			ToolBox(hex).eval_quality(SJ);
			minSJ = 1;
			for (auto c : hex.iter_cells()) minSJ = std::min(SJ[c], minSJ);
			aveSJ = 0;
			for (auto c : hex.iter_cells()) aveSJ += SJ[c]/hex.ncells();
		}

		{// eval max_dist_to_tet
//			
			vec3 closest_pair[2];
			max_dist_to_tet = 0;
			Triangles tri;
			ToolBox(tet).boundary(tri, false);
			tri.connect();
			//DropSurface(tri).apply("tri");
			NearestPointOnTriangles proj_on_tri(tri);
			InsideSurface in_tri(tri);
			for (auto c : hex.iter_cells()) {
				double max_length = 0;
				for (auto f : c.iter_facets()) for (auto h : f.iter_halfedges())
					max_length = std::max(max_length, (h.to().pos() - h.from().pos()).norm2());
				max_length = std::sqrt(max_length);
				int nb_div = std::max(1.0, max_length / eps);
				if (nb_div > 20) {// more than 8000 points to evaluate
					Trace::alert("Too many points to evaluate in an hex");
					nb_div = 20;
				}

				double step = 1. / double(nb_div);

				for (double u = 0; u < 1. + .1 * step; u += step)
					for (double v = 0; v < 1. + .1 * step; v += step)
						for (double w = 0; w < 1. + .1 * step; w += step) {
							vec3 P(0, 0, 0);
							FOR(di, 2)FOR(dj, 2)FOR(dk, 2)
								P += ((di * u) + (1 - di) * (1. - u))
								* ((dj * v) + (1 - dj) * (1. - v))
								* ((dk * w) + (1 - dk) * (1. - w))
								* c.vertex(di + 2 * dj + 4 * dk).pos();
							if (in_tri.request(P)) continue;
							vec3 H = proj_on_tri.request(P);
							if (max_dist_to_tet < (H - P).norm()) {
								closest_pair[0] = P;
								closest_pair[1] = H;
								max_dist_to_tet = (H - P).norm();
							}
						}
			}
			if (show_max_dist  ) drop_arrow(closest_pair[0], closest_pair[1], "max_dist_to_tet");

		}




		tet.connect();
		{
			// eval max_dist_to_hex
			max_dist_to_hex = 0;
			vec3 closest_pair[2];
			Quads quads;
			ToolBox(hex).boundary(quads, false);
			//quads.connect();
			Triangles tri;
			ToolBox(tri.points).copy_from(quads.points);
			int off_v = tri.points.create_points(quads.nfacets());
			tri.create_facets(4 * quads.nfacets());
			for(auto q:quads.iter_facets()) FOR(lh, 4) {
				tri.points[off_v + q] = Quad3(q).bary_verts();
				tri.vert(4 * q + lh, 0) = quads.vert(q, lh);
				tri.vert(4 * q + lh, 1) = quads.vert(q, (lh + 1) % 4);
				tri.vert(4 * q + lh, 2) = off_v + q;
			}
			tri.connect();
			//DropSurface(tri).apply("tri_quad");
			NearestPointOnTriangles proj_on_tri(tri);
			InsideSurface in_tri(tri);
			for (auto c : tet.iter_cells()) {
				double max_length = 0;
				for (auto f : c.iter_facets()) for (auto h : f.iter_halfedges())
					max_length = std::max(max_length, (h.to().pos() - h.from().pos()).norm2());
				max_length = std::sqrt(max_length);
				int nb_div = std::max(1.0, max_length / eps);
				if (nb_div > 20) {// more than 8000 points to evaluate
					Trace::alert("Too many points to evaluate in a tet");
					nb_div = 20;
				}
				double step = 1. / double(nb_div);


				vec3 O = c.vertex(0).pos();
				vec3 b[3];
				FOR(lv, 3) b[lv] = c.vertex(lv + 1).pos() - O;

				for (double u = 0; u < 1. + .1 * step; u += step)
					for (double v = 0; v + u < 1. + .1 * step; v += step)
						for (double w = 0; w + u + v < 1. + .1 * step; w += step) {
							vec3 P = O + u * b[0] + v * b[1] + w * b[2];
							if (in_tri.request(P)) continue;
							vec3 H = proj_on_tri.request(P);
							if (max_dist_to_hex < (H - P).norm()) {
								closest_pair[0] = P;
								closest_pair[1] = H;
								max_dist_to_hex = (H - P).norm();
							}
						}
			}
			if (show_max_dist ) drop_arrow(closest_pair[0], closest_pair[1], "max_dist_to_hex");
		}
		vec3 bbox_size = ToolBox(tet.points).bbox().max-ToolBox(tet.points).bbox().min;
		double max_length=0;
		FOR(d,3) max_length = std::max(max_length ,bbox_size[d]);
		max_dist_to_tet /= max_length;
		max_dist_to_hex /= max_length;
	}
	void HexRemeshQuality::drop_on_trace(std::string prefix ) {
		Trace::log_value(prefix + "_minSJ", minSJ);
		Trace::log_value(prefix + "_aveSJ", aveSJ);
		Trace::log_value(prefix + "_bbox_normalized_dist2tet", max_dist_to_tet);
		Trace::log_value(prefix + "_bbox_normalized_dist2hex", max_dist_to_hex);
		Trace::log_value(prefix + "_bbox_normalized_Hausdorff", std::max(max_dist_to_hex, max_dist_to_tet));

		int nb_edges = 0;
		FOR(i, 10) nb_edges += edge_valence[i] ;
		Trace::log_value(prefix + "_#edges", nb_edges);
		Trace::log_value(prefix + "_#bound_edges", edge_valence[0]);
		Trace::log_value(prefix + "_#singular_edges", nb_edges -edge_valence[0]- edge_valence[4]);


	}
	void HexRemeshQuality::drop_on_shell() {
		plop(minSJ);
		plop(max_dist_to_tet);
		plop(max_dist_to_hex);
	}
