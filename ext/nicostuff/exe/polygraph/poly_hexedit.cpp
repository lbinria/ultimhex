#include <ultimaille/all.h>
#include <framework/trace.h>

#include <framework/nico_framework.h>
#include <framework/benjamin_API.h>
#include "dirty/readonly_mesh_extract_3d.h"
#include <algo/surface/pointset_in_surface.h>


using namespace UM::Linear;

std::string path;
bool run_from_graphite;

void check_hex_validity(Hexahedra& hex, CellFacetAttribute<int>& emb, std::string msg) {
	EdgeGraph eg(hex);
	for (auto e : eg.iter_edges()) if (!e.opposite().active()) {
		Trace::abort(msg + ": hex mesh has a non manifold edge");
	}
	for (auto f : hex.iter_facets()) if (f.on_boundary() && emb[f] < 0) {
		Trace::abort(msg + " : hex mesh a facet without embedding");
	}
}

void save_hex_if_valid(Hexahedra& hex, CellFacetAttribute<int>& emb) {
	check_hex_validity(hex, emb, "hex mesh validity test FAILED ");
	if (run_from_graphite) DropVolume(hex).add(emb, "emb")._just_save_filename(path + "/hex.geogram").apply();
}

void old_smooth() {

	Hexahedra hex;
	auto attr = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
	hex.connect();

	Triangles tri;
	read_by_extension(path + "/tri.geogram", tri);
	tri.connect();

	HexBoundary bound(hex);
	Quads& quad = bound.quad;

	Trace::step("Project on surface");
	{
		PointSet bary;
		bary.create_points(quad.nfacets());
		for (auto f : quad.iter_facets()) bary[f] = Quad3(f).bary_verts();
		PointSetEmbedding  emb(bary);
		ToolBox(emb.tri_emb).copy_from(tri, true);
		emb.tri_emb.connect();
		for (auto f : quad.iter_facets())  emb.set_embedding(f, 2, emb_attr[bound.hex_facet(f)]);
		PointAttribute<vec3> dest(bary);
		FOR(v, bary.size()) dest[v] = bary[v];
		emb.project();
		emb.move_toward(dest);


		{// show boundary constraints
			PointAttribute<vec3> n(bary);
			for (auto f : quad.iter_facets())n[f] = -Triangle3(Surface::Facet(tri, emb.id[f])).normal();;
			Drop(bary, n).apply("n");
			DropPointSet(bary).apply("bary");
			DropPointSet(emb.pts).apply("emb.pts");
		}
		Trace::step("Project/optim");
		LeastSquares ls(3 * quad.nverts());
		// poisson part
		//for (auto h : hex.iter_halfedges()) FOR(d, 3)
		//	ls.add_to_energy(X(d + 3 * h.to()) - X(d + 3 * h.from()) - h.geom().vector()[d]);



		// LAPLACIEN == no long edges
		double eps = 0;
		for (auto h : hex.iter_halfedges()) FOR(d, 3)
			ls.add_to_energy(eps * (X(d + 3 * h.to()) - X(d + 3 * h.from())));

		// h == next h in same direction
		eps = 1;
		for (auto h : hex.iter_halfedges()) {
			auto prev = h.opposite_f().prev();
			auto next = h.opposite_c();
			if (!next.active()) continue;
			next = next.opposite_f().next();
			FOR(d, 3)
				ls.add_to_energy(eps * (
					X(d + 3 * prev.to()) - X(d + 3 * prev.from())
					- X(d + 3 * next.to()) + X(d + 3 * next.from())
					));
		}

		// more or less fix first layer
		{
			double ave = ToolBox(hex).ave_edge_size();
			for (auto f : quad.iter_facets()) {
				Surface::Facet tri_f(tri, emb.id[f]);
				vec3 n = Triangle3(tri_f).normal();
				for (auto tri_h : f.iter_halfedges()) {
					Volume::Halfedge h = bound.hex_halfedge(tri_h).opposite_f().next();
					// the dot
					double w = 5;
					if (true) {//dot
						LinExpr line = .3 * ave;
						FOR(d, 3) line += n[d] * (X(d + 3 * h.from()) - X(d + 3 * h.to()));
						ls.add_to_energy(w * line);
					}
					else {// pos
						FOR(d, 3) ls.add_to_energy(
							w * (X(d + 3 * h.from()) - X(d + 3 * h.to()) + .3 * ave * n[d])
						);
					}
				}
			}
		}


		// parallelograms
		for (auto h : hex.iter_halfedges()) {
			double w = .1;
			FOR(d, 3)
				ls.add_to_energy(
					w * (X(d + 3 * h.to()) - X(d + 3 * h.from())
						+ X(d + 3 * h.next().next().to()) - X(d + 3 * h.next().next().from())
						));
		}

		// fitting
		for (auto f : quad.iter_facets()) {
			Surface::Facet tri_f(tri, emb.id[f]);
			vec3 n = Triangle3(tri_f).normal();
			for (auto h : f.iter_halfedges()) {
				LinExpr line;
				FOR(d, 3) line += n[d] * X(d + 3 * h.from()) - n[d] * bary[f][d];
				ls.add_to_energy(10. * line);
			}
		}
		ls.solve();
		DropVolume(hex).apply("raw hex");
		for (auto v : quad.iter_vertices()) FOR(d, 3) v.pos()[d] = ls.value(d + 3 * v);
		DropVolume(hex).apply("smoothhex");
	}

	save_hex_if_valid(hex, emb_attr);
}

void clear_layer_kill() {
	Hexahedra hex;
	auto attr = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
	hex.connect();
	DropVolume(hex).add(emb_attr, "emb")._just_save_filename(path + "/hex.geogram").apply();
}


void shortest_path_on_vertices(EdgeGraph &eg, int src, int dest,std::vector<int>& result) {
	PointAttribute<double> dist(eg, 1e20);
	PointAttribute<int> prev(eg, -1);
	auto cmp = [&](int i, int j) {// -1 needs special treatment for getting the first eltement
		if (i == -1) return true;
		if (j == -1) return false;
		return dist[i] < dist[j];
		};
	std::set<int, decltype(cmp)> heap(cmp);
	dist[src] = 0;
	heap.insert(src);
	while (true) {
		if (heap.empty()) { Trace::alert("Path does not exist"); return; }
		auto cur = heap.upper_bound(-1);
		PolyLine::Vertex v(eg, *cur);
		heap.erase(cur);

		if (v == dest) {
			if (Trace::trace_steps_active) Trace::alert("path found");
			result.push_back(v);
			while (v != src) {
				v = PolyLine::Edge(eg, prev[v]).from();
				result.push_back(v);
			}
			break;
		}

		for (auto cir : v.iter_edges()) {
			auto opp_v = cir.to();
			double new_dist = dist[v] + (v.pos() + opp_v.pos()).norm();

			if (prev[opp_v] == -1 || dist[opp_v] > new_dist) {
				dist[opp_v] = new_dist;
				prev[opp_v] = cir;
				heap.insert(opp_v);
			}
		}
	}
}


void layer_kill() {
	PolyLine pl;
	read_by_extension(path + "/hex.geogram", pl);
	pl.connect();
	Hexahedra hex;
	auto attr = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
	hex.connect();
	Drop(hex, emb_attr).apply("hex");


	// split polyline edges on hex edges when there is an obvious solution
	{
		EdgeGraph eg(hex);
		DropPolyLine(pl).apply("in_pl");
		std::vector<int> edges2add;
		for (auto e : pl.iter_edges()) {
			PolyLine::Vertex from(eg, e.from());
			PolyLine::Vertex to(eg, e.to());


			std::vector<int> verts;

			shortest_path_on_vertices(eg, from, to, verts);

			//auto v = from;			
			//verts.push_back(v);
			//while (true) {
			//	plop(int(v));
			//	auto next_v = v;
			//	for (auto e : v.iter_edges())
			//		if ((e.to().pos() - to.pos()).norm2() < (next_v.pos() - to.pos()).norm2())
			//			next_v = e.to();
			//	verts.push_back(next_v);
			//	if (next_v == v) break;
			//	if (to == next_v) break;
			//	v = next_v;
			//}


			if (verts.size() > 2)FOR(i, verts.size() - 1) {
				edges2add.push_back(verts[i]);
				edges2add.push_back(verts[i + 1]);
			}
		}
		pl.disconnect();
		int off_e = pl.create_edges(edges2add.size() / 2);
		FOR(le, edges2add.size() / 2) FOR(i, 2) pl.vert(off_e + le, i) = edges2add[2 * le + i];
		pl.connect();

		DropPolyLine(pl).apply("split_pl");
	}

	PointAttribute<int> map(hex, -1);
	{
		EdgeGraph eg(hex);
		EdgeAttribute<int> hex_layer(eg);
		// compute hex layers
		{
			DisjointSet ds(eg.nedges());
			for (auto h : hex.iter_halfedges())
				ds.merge(eg.edge_from_halfedge(h), eg.edge_from_halfedge(h.next().next().opposite_f()));
			int nlayers = ds.get_sets_id(hex_layer.ptr->data);
			Drop(eg, hex_layer).apply("hexlayer");
		}

		EdgeAttribute<bool> collapsable(eg, true);
		// compute edges that can be collapsed
		{
			for (auto h : hex.iter_halfedges()) {
				if (!h.facet().on_boundary()) continue;
				auto h_perp = h.opposite_f().next();
				std::vector<int> in_layer;
				for (auto cir : h_perp.iter_CCW_around_edge()) {
					in_layer.push_back(cir.prev().from());
					in_layer.push_back(cir.prev().opposite_f().prev().from());
				}
				for (auto e : eg.vertex(h.from()).iter_edges()) {
					bool is_in_layer = false;
					for (int v : in_layer)is_in_layer = is_in_layer || (e.to() == v);
					collapsable[e] = collapsable[e] && is_in_layer;
				}

			}
			Drop(eg, collapsable).apply("uncollapsable");
		}

		EdgeAttribute<bool> collapse(eg, false);
		for (auto e_pl : pl.iter_edges()) {
			int layer = -1;

			// find the layer to remove
			{
				PolyLine::Vertex from(eg, e_pl.from());
				for (auto e : from.iter_edges()) {
					if (e.to() != e_pl.to()) continue;
					layer = hex_layer[e];
				}
			}
			if (layer == -1) { Trace::alert("Edge to collapse is not an hex edge"); continue; }

			for (auto h : hex.iter_halfedges()) {
				auto e = eg.edge_from_halfedge(h);
				if (layer != hex_layer[e]) continue;
				if (!collapsable[e] && !collapsable[e.opposite()]) layer = -1;
			}

			if (layer == -1) { Trace::alert("Cannot collapse a self intersecting hex layer"); continue; }

			EdgeAttribute<bool> test_collapse(eg);
			for (auto e : eg.iter_edges()) test_collapse[e] = collapse[e];

			bool can_collapse_layer = true;
			for (auto e : eg.iter_edges()) if (hex_layer[e] == layer) {
				if (collapsable[e]) test_collapse[e] = true;
				else if (collapsable[e.opposite()]) test_collapse[e.opposite()] = true;
				else can_collapse_layer = false;
			}
			if (!can_collapse_layer) { Trace::alert("Cannot collapse the hex layer"); continue; }


			for (auto v : eg.iter_vertices()) {
				std::set<int> dest;
				dest.insert(v);

				bool done = false;
				while (!done) {
					//	plop(v);

					done = true;
					std::vector<int> prev;
					for (auto it : dest) prev.push_back(it);
					dest.clear();
					for (int vid : prev) {
						auto v = eg.vertex(vid);
						bool is_final = true;
						for (auto it : v.iter_edges()) if (test_collapse[it]) {
							dest.insert(it.to());
							is_final = false;

						}
						if (is_final) dest.insert(v);
						else done = false;
					}

				}
				if (dest.size() > 1) can_collapse_layer = false;
			}


			if (!can_collapse_layer) { Trace::alert("A vertex want to move to more than one other vertex"); continue; }

			for (auto e : eg.iter_edges()) collapse[e] = test_collapse[e];
		}


		for (auto v : hex.iter_vertices()) map[v] = v;
		{
			bool done = false;
			while (!done) {
				done = true;
				for (auto e : eg.iter_edges()) if (collapse[e]) if (map[e.from()] != map[e.to()]) {
					done = false;
					map[e.from()] = map[e.to()];
				}
			}
		}
	}

	//WARNING: hex is connected when we modify its verts... not very clean, but quite convinient for embedding propagation
	FOR(c, hex.ncells()) FOR(lv, 8) hex.vert(c, lv) = map[hex.vert(c, lv)];

	std::vector<bool>  to_kill(hex.ncells(), false);
	FOR(c, hex.ncells()) FOR(lv1, 8)FOR(lv2, lv1)
		to_kill[c] = to_kill[c] || (hex.vert(c, lv1) == hex.vert(c, lv2));

	//Drop(hex, emb_attr).apply("emb"); return;
	//for (auto f : hex.iter_facets()) if (!f.on_boundary())emb_attr[f] = -1;
	{// propagate embedding
		bool done = false;
		while (!done) {
			done = true;
			for (auto f : hex.iter_facets()) {
				if (emb_attr[f] < 0) continue;
				if (!to_kill[f.cell()]) continue;
				auto f_next = f.halfedge(0).opposite_f().next().next().opposite_f().facet().opposite();
				if (f_next.active()) {
					std::swap(emb_attr[f], emb_attr[f_next]);
					done = false;
				}
			}
		}
	}
	hex.disconnect();


	hex.delete_cells(to_kill);
	hex.delete_isolated_vertices();
	hex.connect();

	for (auto f : hex.iter_facets()) if (f.on_boundary() && emb_attr[f] < 0) {
		for (auto in_f : hex.iter_facets()) if (emb_attr[in_f] >= 0 || !in_f.on_boundary()) emb_attr[in_f] = 0;
		Drop(hex, emb_attr)._skip_value(0).apply("feeedback");
		Drop(hex, emb_attr)._skip_value(0)._just_save_filename(path + "/feedback.geogram").apply();
		return;
	}
	Drop(hex, emb_attr).apply("hex");
	save_hex_if_valid(hex, emb_attr);
}

void pad_interface() {
	Hexahedra hex;
	read_by_extension(path + "/hex.geogram", hex);
	hex.connect();

	PointAttribute<bool> on_border(hex, false);
	for (auto h : hex.iter_halfedges()) if (h.facet().on_boundary()) on_border[h.from()] = true;

	{
		Quads quad;
		quad.points = hex.points;
		for (auto f : hex.iter_facets()) {
			if (!f.on_boundary())  continue;
			int qf = quad.create_facets(1);
			FOR(lv, 4) quad.vert(qf, lv) = f.vertex(lv);
		}
		FacetAttribute<int> fpaint(quad, 0);
		DropSurface(quad).add(fpaint, "fpaint")._just_save_filename(path + "/padselecton.geogram").apply();
		Drop(quad, fpaint).apply("pad_interface");
	}
	{
		Quads quad;
		quad.points = hex.points;
		for (auto f : hex.iter_facets()) {
			if (f.on_boundary())  continue;
			// how to continue the outerloop inside a inner loop ?
			if ([&]() {FOR(lv, 4)  if (on_border[f.vertex(lv)]) return false;  return true;  }()) continue;
			int qf = quad.create_facets(1);
			FOR(lv, 4) quad.vert(qf, lv) = f.vertex(lv);
		}
		FacetAttribute<int> fpaint(quad, 0);
		DropSurface(quad).add(fpaint, "fpaint")._just_save_filename(path + "/padselectin.geogram").apply();

		Drop(quad, fpaint).apply("pad_interface");
	}
}

void compute_padding_faces(Hexahedra& hex, CellFacetAttribute<bool>& pad_face,bool invert_faces) {
	CellFacetAttribute<int> constraint(hex, 0);
	hex.connect();

	for (auto filename : { "/padselectin.geogram" ,"/padselecton.geogram" }) { // load attribute from quads
		// TODO: replace quads by another surface with some inner quads of the hex mesh => for inner constraints
		Quads load_quad;
		auto attr = read_by_extension(path + filename, load_quad);
		load_quad.connect();
		FacetAttribute<int> load_fpaint("fpaint", attr, load_quad);
		Drop(load_quad, load_fpaint).apply("inon");
		for (auto f : hex.iter_facets()) {
			Surface::Vertex v(load_quad, f.vertex(0));
			for (auto h : v.iter_halfedges()) {
				if (h.to() == f.vertex(1)
					&& h.next().to() == f.vertex(2)
					&& h.next().next().to() == f.vertex(3)
					)
					if (load_fpaint[h.facet()] == 1) constraint[f] = 1;
			}
		}
	}


	Drop(hex, constraint).apply("input constraints");

	CellFacetAttribute<int> layer(hex);
	{
		DisjointSet ds(hex.nfacets());
		for (auto h : hex.iter_halfedges()) {
			auto opp = h.opposite_c();
			if (!opp.active()) continue;
			ds.merge(h.opposite_f().facet(), opp.opposite_f().facet());
		}
		ds.get_sets_id(layer.ptr->data);
	}
	Drop(hex, layer).apply("layers");

	std::set<int> active_layers;
	for (auto f : hex.iter_facets()) if (constraint[f] == 1) active_layers.insert(layer[f]);
	for (auto l : active_layers) {
		for (auto f : hex.iter_facets()) if (layer[f] == l && constraint[f] != 1) constraint[f] = 2;
	}
	Drop(hex, constraint)._skip_value(0).apply("input constraints");


	EdgeGraph eg(hex);
	EdgeAttribute<int> ncharts(eg, 0);
	for (auto h : hex.iter_halfedges())
		if (constraint[h.facet()] > 0)
			ncharts[eg.edge_from_halfedge(h)]++;
	for (auto e : eg.iter_edges()) if (e > e.opposite()) {
		ncharts[e] += ncharts[e.opposite()];
		ncharts[e.opposite()] = ncharts[e];
	}
	Drop(eg, ncharts)._skip_value(0).apply_wireframe("cuts");
	if (run_from_graphite)Drop(eg, ncharts)._skip_value(0)._just_save_filename(path + "/cuts.geogram").apply_wireframe("cuts");
	CellFacetAttribute<int> layer_cut(hex);
	{
		DisjointSet ds(hex.nfacets());
		for (auto h : hex.iter_halfedges()) {
			auto opp = h.opposite_c();
			if (!opp.active()) continue;
			if (ncharts[eg.edge_from_halfedge(h)] > 2) continue;
			ds.merge(h.opposite_f().facet(), opp.opposite_f().facet());
		}
		ds.get_sets_id(layer_cut.ptr->data);
	}
	Drop(hex, layer_cut).apply("layer_cut");


	CellFacetAttribute<int> constraint_cut(hex,-1);

	std::set<int> active_layers_cut;
	for (auto f : hex.iter_facets()) if (constraint[f] == 1) active_layers_cut.insert(layer_cut[f]);
	for (auto l : active_layers_cut) {
		for (auto f : hex.iter_facets()) if (layer_cut[f] == l) constraint_cut[f] = l;
	}
	Drop(hex, constraint_cut)._skip_value(-1).apply("constraints cut");

	for (auto f : hex.iter_facets()) pad_face[f] = constraint_cut[f] >= 0;


	if (invert_faces) {
		CellFacetAttribute<bool> tmp(hex);
		for (auto f : hex.iter_facets()) tmp[f] = pad_face[f];
		for (auto f : hex.iter_facets()) pad_face[f]=false;
		for (auto f : hex.iter_facets()) if (tmp[f]) {
			//pad_face[f] = false;
			if (f.opposite().active()) pad_face[f.opposite()] = true;
		}
	}
}

void pad_faces(bool invert_faces) {
	Hexahedra hex;
	auto attr = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
	hex.connect();
	check_hex_validity(hex, emb_attr, "CHECK input hex mesh");

	CellFacetAttribute<bool> pad_face(hex);
	compute_padding_faces(hex, pad_face, invert_faces);

	EdgeGraph eg(hex);
	EdgeAttribute<int> nfaces(eg, 0);
	for (auto h : hex.iter_halfedges()) if (pad_face[h.facet()]) nfaces[eg.edge_from_halfedge(h)]++;
	//Drop(eg, nfaces)._skip_value(0).apply_wireframe("nfaces");
	PolyLine pl;
	pl.points = hex.points;
	for (auto e : eg.iter_edges()) if (nfaces[e] > 1) ToolBox(pl).add_segment(e.from(),e.to());
	


	EdgeAttribute<bool> boundary(eg, false);
	EdgeAttribute<int> valence(eg, 0);
	for (auto h : hex.iter_halfedges()) if (h.facet().on_boundary()) boundary[eg.edge_from_halfedge(h)] = true;

	for (auto e : eg.iter_edges()) if (!boundary[e] && nfaces[e] != nfaces[e.opposite()] ) ToolBox(pl).add_segment(e.from(), e.to());


	for (auto h : hex.iter_halfedges()) valence[eg.edge_from_halfedge(h)]++;
	for (auto h : hex.iter_halfedges()) 
		if (
			h.facet().on_boundary() 
			&& pad_face[h.facet()]
			&&  nfaces[eg.edge_from_halfedge(h)] == 1 
			&& nfaces[eg.edge_from_halfedge(h).opposite()] != 1
			&& valence[eg.edge_from_halfedge(h)]>1) {
		ToolBox(pl).add_segment(h.from(), h.to());
	}
	
	DropPolyLine(pl).apply();
	if (run_from_graphite)DropPolyLine(pl)._just_save_filename(path + "/feedbackerror.geogram").apply();

	Drop(hex, pad_face)._skip_value(false).apply("feeedback"); 
	if (run_from_graphite)Drop(hex, pad_face)._skip_value(false)._just_save_filename(path + "/feedback.geogram").apply();
}

void pad_from_fpaint(int nb_layers , bool invert_faces) {
	Hexahedra hex;
	auto attr = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
	hex.connect();
	check_hex_validity(hex, emb_attr, "CHECK input hex mesh");

	CellFacetAttribute<bool> pad_face(hex);
	compute_padding_faces(hex, pad_face, invert_faces);

	//BenjaminAPI::pad(hex, pad_face);
	FOR(it, nb_layers) {
		HexPad padder(hex);
		padder.apply(pad_face);
		DropVolume(hex).apply("padded");
		{// upate emb

			// update boundary of inside layer
			for (auto h : hex.iter_halfedges()) {
				if (!pad_face[h.facet()]) continue;
				auto opp = h.opposite_c();
				um_assert(opp.active()); // need to check how it can happen after the padding...
				auto new_face = opp.opposite_f().facet();
				if (!new_face.on_boundary()) continue;
				emb_attr[new_face] = emb_attr[h.opposite_f().facet()];
			}
			// update boundary layer
			for (auto h : hex.iter_halfedges()) {
				if (!pad_face[h.facet()]) continue;
				auto opp = h.opposite_c();
				um_assert(opp.active()); // need to check how it can happen after the padding...
				auto new_face = opp.opposite_f().next().next().opposite_f().facet();
				if (!new_face.on_boundary()) continue;
				emb_attr[new_face] = emb_attr[h.facet()];
			}
		}
	}
	// clear emb inside the volume
	for (auto f : hex.iter_facets()) if (!f.on_boundary())emb_attr[f] = -1;

	Drop(hex, emb_attr).apply("emb_attr");
	Drop(hex, pad_face)._skip_value(false).apply("pad faces");
	save_hex_if_valid(hex, emb_attr);
}

void hexview() {
	Hexahedra hex;
	read_by_extension(path + "/hex.geogram", hex);
	hex.connect();
	EdgeGraph eg(hex);
	EdgeAttribute<int> singu(eg, 0);
	for (auto h : hex.iter_halfedges()) singu[eg.edge_from_halfedge(h)] += h.facet().on_boundary() ? 3 : 1;

	Drop(eg, singu)._skip_value(4).apply("singu");

	CellAttribute<int> block(hex);
	CellFacetAttribute<int> layer(hex);
	int nlayers;
	{
		DisjointSet ds(hex.nfacets());
		for (auto h : hex.iter_halfedges()) {
			auto opp = h.opposite_c();
			if (!opp.active()) continue;
			ds.merge(h.opposite_f().facet(), opp.opposite_f().facet());
		}
		nlayers = ds.get_sets_id(layer.ptr->data);
		Drop(hex, layer).apply("layers");
	}

	std::vector<bool> cut(nlayers, false);
	for (auto h : hex.iter_halfedges()) {
		if (singu[eg.edge_from_halfedge(h)] != 4 && !h.facet().on_boundary())
			cut[layer[h.facet()]] = true;
	}
	DisjointSet ds(hex.ncells());
	for (auto f : hex.iter_facets()) {
		auto opp = f.opposite();
		if (!opp.active()) continue;
		if (cut[layer[f]] || cut[layer[opp]]) continue;
		ds.merge(f.cell(), opp.cell());
	}
	ds.get_sets_id(block.ptr->data);
	Drop(hex, block).apply("block");

	if (run_from_graphite) {
		PolyLine pl;
		EdgeAttribute<int> plsingu(pl);
		pl.points = eg.points;
		for (auto e : eg.iter_edges()) if (singu[e] != 4) {
			auto ne = pl.create_edges(1);
			pl.vert(ne, 0) = e.from();
			pl.vert(ne, 1) = e.to();
			plsingu[ne] = singu[e];
		}
		//Drop(pl, plsingu).apply("plsingu");
		DropPolyLine(pl).add(plsingu, "singu")._just_save_filename(path + "/singuview.geogram").apply();
		DropVolume(hex).add(block, "block")._just_save_filename(path + "/hexview.geogram").apply();
	}
}


void embeditinit(bool gmsh_chart) {
	Triangles tri;
	auto tri_attr = read_by_extension(path + "/tri.geogram", tri);
	tri.connect();

	Hexahedra hex;
	auto attr = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr, hex, -1);
	hex.connect();

	FacetAttribute<int> tri_chart(tri);
	if (gmsh_chart) {
		FacetAttribute<int> gmsh_tri_chart("region", tri_attr, tri);
		for (auto f : tri.iter_facets()) tri_chart[f] = gmsh_tri_chart[f];
	}
	else
	{
		CornerAttribute<int> feature(tri);
		FeatureEdgeDetector(tri).dihedral_angle().threshold().remove_small_features().remove_small_features().remove_small_features().apply(feature, false);
		Drop(tri, feature)._wireframe(true).apply_half_edge("features");

		DisjointSet ds(tri.nfacets());
		for (auto h : tri.iter_halfedges()) {
			auto opp = h.opposite();
			um_assert(opp.active());
			if (feature[h] == -1)
				ds.merge(h.facet(), opp.facet());
		}
		ds.get_sets_id(tri_chart.ptr->data);
	}
	Drop(tri, tri_chart).apply("tri_chart");

	HexBoundary bound(hex);
	Quads& quad = bound.quad;
	FacetAttribute<int> quad_chart(quad);
	for (auto f : quad.iter_facets()) 
		quad_chart[f] = tri_chart[emb_attr[bound.hex_facet(f)]];
	
	Drop(quad, quad_chart).apply("tri_chart");
	if (run_from_graphite) {
		DropSurface(tri).add(tri_chart, "chart")._just_save_filename(path + "/trichart.geogram").apply();
		DropSurface(quad).add(quad_chart, "chart")._just_save_filename(path + "/quadchart.geogram").apply();
	}
}


void embeditapply() {
	Hexahedra hex;
	auto attr_hex = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr_hex, hex, -1);
	hex.connect();

	Quads quad;
	auto attr_quad = read_by_extension(path + "/quadchart.geogram", quad);
	FacetAttribute<int> quad_chart("chart", attr_quad, quad, -1);
	quad.connect();

	Triangles tri;
	auto attr_tri= read_by_extension(path + "/trichart.geogram", tri);
	FacetAttribute<int> tri_chart("chart", attr_tri, tri, -1);
	tri.connect();

	Drop(quad, quad_chart).apply("quadchart");
	Drop(tri, tri_chart).apply("trichart");

	FacetAttribute<int> quad2hex_face(quad,-1);
	{
		for (auto f : hex.iter_facets()) {
			if (!f.on_boundary()) continue;
			Surface::Vertex v(quad, f.halfedge(0).from());
			for (auto h : v.iter_halfedges()) if (h.to() == f.halfedge(0).to())
				quad2hex_face[h.facet()] = f;
		}
		for (auto f : quad.iter_facets()) um_assert(quad2hex_face[f] != -1);
	}

	FacetAttribute<int> quad2tri_face(quad, -1);
	FacetAttribute<vec3> bary_pos(quad);
	{
		PolyLine pl;
		for (auto f_quad : quad.iter_facets()) {
			int chart = quad_chart[f_quad];
			vec3 G= Quad3(f_quad).bary_verts();
			quad2tri_face[f_quad] = -1;
			double best_dist2 = 1e20;
			for (auto f_tri : tri.iter_facets()) {
				if (tri_chart[f_tri] != chart) continue;
				vec3 bc = Triangle3(f_tri).bary_coords(G);
				FOR(lv, 3) bc[lv] = std::max(.001, std::min(.999, bc[lv]));
				double sum = bc[0] + bc[1] + bc[2];
				FOR(lv, 3) bc[lv] /= sum;
				vec3 proj = bc[0] * f_tri.vertex(0).pos() + bc[1] * f_tri.vertex(1).pos() + bc[2] * f_tri.vertex(2).pos();
				if ((proj - G).norm2() < best_dist2) {
					quad2tri_face[f_quad] = f_tri;
					bary_pos[f_quad] = proj;
					best_dist2 = (proj - G).norm2();
				}
			}
			ToolBox(pl).add_segment(G, bary_pos[f_quad]);
			um_assert(quad2tri_face[f_quad] != -1);
			emb_attr[quad2hex_face[f_quad]] = quad2tri_face[f_quad];
		}
		DropPolyLine(pl).apply("match");
	}


	save_hex_if_valid(hex, emb_attr);
	return;


	Drop(quad, quad_chart).apply("quadchart");
	enum {BARY,DIFFUSION,LSCM
	} proj_strat = LSCM;
	switch (proj_strat) {
	case BARY:
		for (auto v : quad.iter_vertices()) {
			vec3 sum_P(0, 0, 0);
			double n = 0;
			for (auto h : v.iter_halfedges()) {
				sum_P += bary_pos[h.facet()];
				n += 1;
			}
			v.pos() = sum_P / n;
		}
		break;
	case (DIFFUSION):
		FOR(dim, 3) {
		double t = 10;
		LeastSquares ls(quad.nverts());
			for (auto v : quad.iter_vertices()) {
				LinExpr line;
				if (!v.halfedge().active()) continue;
				for (auto h : v.iter_halfedges())
					line += X(h.to()) - X(v);
				ls.add_to_energy((X(v) - t * line) - v.pos()[dim]);
			}
			ls.solve();
			for (auto v : quad.iter_vertices()) v.pos()[dim] = ls.value(v);
		}
		break;
	case (LSCM):
		{
			LeastSquares ls(3*quad.nverts());
			for (auto f_quad : quad.iter_facets()) {
				Surface::Facet f_tri(tri, quad2tri_face[f_quad]);
				vec3 G = Triangle3(f_tri).bary_verts();
				vec3 n = Triangle3(f_tri).normal();
				Quaternion q;
				q.v = n * std::sin(M_PI / 4.);
				q.w = std::cos(M_PI / 4.);
				mat3x3 R = q.rotation_matrix();
				for (auto h : f_quad.iter_halfedges()) {
					//lscm
					auto other = h.prev().opposite();
					FOR(dim, 3) {
						LinExpr line = X(dim + 3 * h.to()) - X(dim + 3 * h.from());
						FOR(d, 3) line += R[dim][d] * (X(d + 3 * other.to()) - X(d + 3 * other.from()));
						ls.add_to_energy(line);
					}

					{
						LinExpr line  = -n *G;
						FOR(d, 3) line += n[d] * X(d + 3 * h.from());
						ls.add_to_energy(10.*line);
					}
				}
				//FOR(lv,4) FOR(d,3) ls.add_to_energy(X(d+3*f_quad.vertex(lv))-G[d]);
			}
			ls.solve();
			for (auto v : quad.iter_vertices()) FOR(d,3)v.pos()[d] = ls.value(d+3*v);
		}
		break;

	}
	Drop(quad, quad_chart).apply("quadchartsmooth");


}


struct EmbeddedHexSmoother {
	EmbeddedHexSmoother(
		Hexahedra& hex,
		CellFacetAttribute<int>& emb_attr,
		Triangles& tri,
		FacetAttribute<int>& tri_chart) : hex(hex), emb_attr(emb_attr), tri(tri), tri_chart(tri_chart),
		bound(hex), quad(bound.quad), emb(bary) {
		bary.create_points(quad.nfacets());
		for (auto f : quad.iter_facets()) bary[f] = Quad3(f).bary_verts();
		CornerAttribute<bool> feature(tri, false);
		for (auto h : tri.iter_halfedges()) feature[h] = (tri_chart[h.facet()] != tri_chart[h.opposite().facet()]);
		emb.init_from_triangles(tri, &feature);
		for (auto f : quad.iter_facets())  emb.set_embedding(f, 2, emb_attr[bound.hex_facet(f)]);

	}

	Hexahedra& hex;
	CellFacetAttribute<int>& emb_attr;

	Triangles& tri;
	FacetAttribute<int>& tri_chart;
	HexBoundary bound;
	Quads &quad;
	PointSet bary;
	PointSetEmbedding  emb;


	void bary2verts(std::string name) {
		for (auto v : quad.iter_vertices()) {
			vec3 sum_P(0, 0, 0);
			double n = 0;
			for (auto h : v.iter_halfedges()) {
				sum_P += bary[h.facet()];
				n += 1;
			}
			if (n > 0) v.pos() = sum_P / n;
		}
		DropSurface(quad).apply(name);
	};

	void show_boundary_constraints(){// show boundary constraints
		PointAttribute<vec3> n(bary);
		for (auto f : quad.iter_facets())n[f] = Triangle3(Surface::Facet(tri, emb.id[f])).normal();;
		Drop(bary, n).apply("n");
		DropPointSet(bary).apply("bary");
		DropPointSet(emb.pts).apply("emb.pts");
	}

	void lscm() {
		LeastSquares ls(3 * quad.nverts());
		for (auto f_quad : quad.iter_facets()) {
			Surface::Facet f_tri(tri, emb.id[f_quad]);
			vec3 G = Triangle3(f_tri).bary_verts();
			vec3 n = Triangle3(f_tri).normal();
			Quaternion q;
			q.v = n * std::sin(M_PI / 4.);
			q.w = std::cos(M_PI / 4.);
			mat3x3 R = q.rotation_matrix();
			for (auto h : f_quad.iter_halfedges()) {
				auto other = h.prev().opposite();
				FOR(dim, 3) {
					LinExpr line = X(dim + 3 * h.to()) - X(dim + 3 * h.from());
					FOR(d, 3) line += R[dim][d] * (X(d + 3 * other.to()) - X(d + 3 * other.from()));
					ls.add_to_energy(line);
				}
				{
					LinExpr line = -n * G;
					FOR(d, 3) line += n[d] * X(d + 3 * h.from());
					ls.add_to_energy(10. * line);
				}
			}
		}
		ls.solve();
		for (auto v : quad.iter_vertices()) FOR(d, 3)v.pos()[d] = ls.value(d + 3 * v);
	}


	void smooth_inside() {
		LeastSquares ls(3 * quad.nverts());
		// fitting
		for (auto f : quad.iter_facets()) for (auto h : f.iter_halfedges()) 
			FOR(d,3) ls.fix(d + 3 * h.from(),h.from().pos()[d]);
		// LAPLACIEN == no long edges
		double eps = 1;
		for (auto h : hex.iter_halfedges()) FOR(d, 3)
			ls.add_to_energy(eps * (X(d + 3 * h.to()) - X(d + 3 * h.from())));
		ls.solve();
		for (auto v : quad.iter_vertices()) FOR(d, 3) v.pos()[d] = ls.value(d + 3 * v);
		DropVolume(hex).apply("smoothhex");
	}

	void apply() {
		
		
		FOR(it, 3) {
			emb.project();
			PointAttribute<vec3> dest(bary);
			for (auto f : quad.iter_facets()) dest[f] = Quad3(f).bary_verts();
			emb.move_toward(dest);
			show_boundary_constraints();
			lscm();
			DropSurface(quad).apply("lscm");
		}

		smooth_inside();
		DropVolume(hex).apply("smooth_inside");

		//PointAttribute<vec3> dest(bary);
		//FOR(v, bary.size()) dest[v] = bary[v] + vec3(0, 3, 0);
		//bary2verts("quad_init");
		//emb.project();
		//bary2verts("quad_proj");
		//emb.move_toward(dest);
		//bary2verts("quad_walk");
		//emb.project();
		//bary2verts("quad_reproj");
		//show_boundary_constraints();
	}

};

void smooth() {

	Hexahedra hex;
	auto attr_hex = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb_attr("emb", attr_hex, hex, -1);
	hex.connect();

	Triangles tri;
	auto attr_tri = read_by_extension(path + "/trichart.geogram", tri);
	FacetAttribute<int> tri_chart("chart", attr_tri, tri, -1);
	tri.connect();

	Drop(tri, tri_chart).apply("trichart");
	DropVolume(hex).apply("in_hex");
	HexBoundary bound(hex);

	EmbeddedHexSmoother smoother(hex, emb_attr, tri, tri_chart);
	smoother.apply();

	save_hex_if_valid(hex, emb_attr);
}



void killhexinit() {
	Hexahedra hex;
	auto attr_hex = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb("emb", attr_hex, hex, -1);
	CellAttribute<int> tokill(hex, 0);
	if (run_from_graphite)
		DropVolume(hex).add(tokill,"tokill").add(emb, "emb")._just_save_filename(path + "/hex.geogram").apply();
}
void killhexapply() {
	Hexahedra hex;
	auto attr_hex = read_by_extension(path + "/hex.geogram", hex);
	CellFacetAttribute<int> emb("emb", attr_hex, hex, -1);
	CellAttribute<int> to_kill("tokill", attr_hex, hex, 0);
	//Drop(hex, emb).apply("emb");
	hex.connect();

	FOR(iter, 10) {
		CellFacetAttribute<int> nvemb(hex, -1);
		// copy emb where it is fine
		for (auto f : hex.iter_facets()) nvemb[f] = emb[f];
		
		// propagated emb from killed quads
		for (auto f : hex.iter_facets()) {
			if (f.on_boundary()) continue;
			//if (!to_kill[f.cell()]) continue;
			if (nvemb[f.opposite()] != -1) continue;
			FOR(le, 4)
				if (emb[f.halfedge(le).opposite_f().facet()] != -1)
					nvemb[f.opposite()] = emb[f.halfedge(le).opposite_f().facet()];
		}

		for (auto f : hex.iter_facets()) emb[f] = nvemb[f];
	}

	Drop(hex, to_kill).apply("tokill");

	std::vector<bool> killvect(hex.ncells());
	FOR(c, hex.ncells()) killvect[c] = (to_kill[c]==1);
	hex.disconnect();
	hex.delete_cells(killvect);
	hex.connect();

	for (auto f : hex.iter_facets()) 
		if (!f.on_boundary()) emb[f]  = -1;
	
	Drop(hex, emb).apply("emb");
	save_hex_if_valid(hex,emb);
}



void framework_parameters(NicoFramework& fw) {

	fw.add("enum", "algo", "embeditinit").possible_values("smooth,clearlayerkill,layerkill,pad,padInterface,padfaces,hexview,embeditinit,embeditapply,killhexapply,killhexinit");
	//fw.add("string", "projectpath", "C:/NICO/data/polycubeprojects/CAD/easy_pont");
	fw.add("string", "projectpath", "C:/NICO/data/mambo-master/mesh/B18.step");
	fw.add("bool", "gmshchart", "true");
	fw.add("bool", "invertNormal", "false");
	
	fw.add("int", "nbpadlayers", "1");

}


void framework_main(NicoFramework& fw) {
	path = fw["projectpath"];
	run_from_graphite = fw["run_from"].is("graphite");

	if (fw["algo"].is("clearlayerkill"))	clear_layer_kill();
	if (fw["algo"].is("layerkill"))			layer_kill();
	if (fw["algo"].is("smooth"))			smooth();

	if (fw["algo"].is("padInterface")) 		pad_interface();
	if (fw["algo"].is("padfaces")) 			pad_faces(fw["invertNormal"]);
	if (fw["algo"].is("pad")) 				pad_from_fpaint(fw["nbpadlayers"], fw["invertNormal"]);

	if (fw["algo"].is("hexview")) 			hexview();
	if (fw["algo"].is("embeditinit")) 		embeditinit(fw["gmshchart"]);
	if (fw["algo"].is("embeditapply")) 		embeditapply();

	if (fw["algo"].is("killhexinit")) 		killhexinit();
	if (fw["algo"].is("killhexapply")) 		killhexapply();




}
