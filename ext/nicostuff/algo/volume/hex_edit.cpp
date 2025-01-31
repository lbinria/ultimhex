
#include <ultimaille/all.h>
#include <framework/trace.h>


#include "toolbox.h"
#include "drop_attribute.h"
#include <algo/volume/hex_select.h>
#include <algo/volume/hex_edit.h>
#include <algo/volume/hex_smooth.h>



HexPad::HexPad(Hexahedra& hex) : hex(hex) {}
void HexPad::apply(CellFacetAttribute<bool>& pad_face, bool verbose, bool geom_optim,bool dilate ,PointAttribute<vec3>* other_pos ) {
	
	if (!hex.connected()) {
		hex.connect();
		plop("HexPad::apply called on not connected mesh");
	}
	if (verbose) Trace::step("group corners");

	CellCornerAttribute<int> next_in_group(hex);
	{
		FOR(c, hex.ncorners())next_in_group[c] = c;

		for (auto h : hex.iter_halfedges()) {
			auto opp = h.opposite_c();
			if (!opp.active()) continue;
			if (pad_face[h.facet()]) continue;
			if (pad_face[opp.facet()]) continue;
			int c0 = h.from_corner();
			int c1 = opp.next().from_corner();
			int it = c0;
			do {
				if (it == c1) break;
				it = next_in_group[it];
			} while (it != c0);
			if (it == c1) continue;
			std::swap(next_in_group[c0], next_in_group[c1]);
		}
	}

	if (verbose) Trace::step("c2v");

	CellCornerAttribute<int> c2v(hex, -1);


	for (auto h : hex.iter_halfedges()) {
		if (!pad_face[h.facet()]) continue; // affected to an original vertex
		int c = h.from_corner();
		if (c2v[c] != -1) continue;			// already done
		int v = hex.points.create_points(1);
		hex.points[v] = h.from().pos();
		if (other_pos!=NULL) (*other_pos)[v] =  (*other_pos)[h.from()];
		int it = c;
		do {
			c2v[it] = v;
			it = next_in_group[it];
		} while (it != c);
	}

	for (auto h : hex.iter_halfedges())
		if (c2v[h.from_corner()] == -1) c2v[h.from_corner()] = h.from();


	for (auto h : hex.iter_halfedges()) {
		if (!pad_face[h.facet()]) continue; // affected to an original vertex
		int c = h.from_corner();
		int v = c2v[c];
		hex.points[v] += .1 * (h.opposite_f().next().to().pos() - h.from().pos());
		//if (other_pos!=NULL) (*other_pos)[v] += .1 * ((*other_pos)[h.opposite_f().next().to()] - (*other_pos)[h.from()]);
	}



	std::vector<std::array<int, 8> > new_hex;
	for (auto f : hex.iter_facets()) {
		if (!pad_face[f]) continue;
		int h[4] = { f.halfedge(0), f.halfedge(1), f.halfedge(3), f.halfedge(2) };
		std::array<int, 8> nh;
		FOR(i, 4) nh[i + 4] = Volume::Halfedge(hex, h[i]).from();
		FOR(i, 4) nh[i] = c2v[Volume::Halfedge(hex, h[i]).from_corner()];
		new_hex.push_back(nh);
	}



	if (verbose) Trace::step("rebind old to vertices");
	FOR(c, hex.ncorners()) hex.vert(c / 8, c % 8) = c2v[c];

	if (verbose) Trace::step("create new hex");

	int off_c = hex.create_cells(new_hex.size());
	FOR(i, new_hex.size()) FOR(lv, 8) hex.vert(off_c + i, lv) = new_hex[i][lv];
	if (verbose) Drop(hex, c2v).apply("geom pad");

	hex.connect();
	if (geom_optim) {
		plop("Coucou benjamin");
		double ave = ToolBox(hex).ave_edge_size();
		for (auto v : hex.iter_vertices()) v.pos()/=ave;

		HexSmoother smoother(hex);

		HexPointSelect select(hex);
		select.fill(true);
		CellFacetAttribute<int> chart(hex,0);
		for (auto f : hex.iter_facets()) if (pad_face[f]) chart[f] = 1;
		select.by_facet_value(chart, 1, false);
		if (dilate) select.dilate(false);		

		for (auto f : hex.iter_facets()) {
			auto opp = f.opposite();
			if (opp.active()) continue;
			vec3 n = ToolBox(hex).facet_geom(f).normal();
			FOR(lv, 4) smoother.add_constraint(f.vertex(lv),f.vertex(lv).pos(),n );
		}

		smoother.set_lock(select);
		smoother.smooth_elliptic();

		for (auto v : hex.iter_vertices()) v.pos() *= ave;
	}


}



void refine_marked_cells(Hexahedra& hex, CellAttribute<bool>& need_refine,bool geom_optim ) {
	hex.connect();
	// mark halfedges
	std::vector<bool> halfedge2refine(hex.nfacets() * 4, false);

	auto check_halfedge2refine = [&]() {
		for (auto h : hex.iter_halfedges()) if (h.opposite_c().active())
			um_assert(halfedge2refine[h] == halfedge2refine[h.opposite_c()]);
		for (auto h : hex.iter_halfedges())
			um_assert(halfedge2refine[h] == halfedge2refine[h.opposite_f()]);
	};


	for (auto h : hex.iter_halfedges()) if (need_refine[h.facet().cell()])
		for (auto cir : h.iter_CCW_around_edge())
			halfedge2refine[cir] = true;
	check_halfedge2refine();

	auto show_edge2refine = [&]() {
		PolyLine pl;
		pl.points = hex.points;
		for (auto h : hex.iter_halfedges())
			if (halfedge2refine[h]) {
				int s = pl.create_edges(1);
				pl.vert(s, 0) = h.from();
				pl.vert(s, 1) = h.to();
			}
		DropPolyLine(pl).apply("pl");
	};

	for (auto seed_h : hex.iter_halfedges()) {
		if (seed_h >= halfedge2refine.size()) break;
		if (!halfedge2refine[seed_h]) continue;
		std::vector<bool> h_in_hex_layer(hex.nfacets() * 4, false);


		{
			std::vector<int> layer(hex.nfacets() * 4, -1);
			int nb_layers = 0;
			for (auto seed_layer : hex.iter_halfedges()) {
				if (seed_layer >= halfedge2refine.size()) break;
				if (!halfedge2refine[seed_layer]) continue;
				if (layer[seed_layer] != -1) continue;
				std::vector<int> stack;
				layer[seed_layer] = nb_layers;
				stack.push_back(seed_layer);
				while (!stack.empty()) {
					Volume::Halfedge h(hex, stack.back());
					stack.pop_back();

					// turn around
					for (auto cir : h.iter_CCW_around_edge()) {
						if (layer[cir] != -1) continue;
						layer[cir] = nb_layers;
						stack.push_back(cir);
					}
					// go next
					auto next = h.next().next().opposite_f();
					if (halfedge2refine[next] && layer[next] != nb_layers) {
						layer[next] = nb_layers;
						stack.push_back(next);
					}
				}
				nb_layers++;
			}

			std::vector<bool> active_layers(nb_layers, false);
			active_layers[0] = true;
			FOR(l, nb_layers) {
				active_layers[l] = true;
				for (auto h : hex.iter_halfedges()) {
					if (layer[h] != l) continue;
					if (layer[h.opposite_f()] < layer[h]) active_layers[l] = false;
					if (layer[h.next()] != -1)
						if (active_layers[layer[h.next()]]) active_layers[l] = false;

					if (layer[h.prev()] != -1)
						if (active_layers[layer[h.prev()]]) active_layers[l] = false;
				}
			}
			for (auto h : hex.iter_halfedges()) {
				if (layer[h] == -1) continue;
				h_in_hex_layer[h] = active_layers[layer[h]];
			}
		}



		// mark pad faces from halfedges to be refined
		CellFacetAttribute<bool> pad_face(hex, false);
		{
			CellFacetAttribute<int> pad_face_type(hex, -1);
			for (auto h : hex.iter_halfedges()) {
				if (!h_in_hex_layer[h]) continue;
				pad_face_type[h.prev().opposite_f().facet()] = 0;
				pad_face_type[h.next().opposite_f().facet()] = 0;
				auto cir = h.prev().opposite_f();
				do {
					if (!h_in_hex_layer[cir.opposite_f().prev().opposite_f()]
						&& !h_in_hex_layer[cir.opposite_f().next()]
						&& cir.opposite_f().opposite_c().active()
						)
						pad_face_type[cir.opposite_f().facet()] = 1;
					cir = cir.next();
				} while (cir != h.prev().opposite_f());
			}
			for (auto f : hex.iter_facets()) {
				auto opp = f.opposite();
				if (!opp.active()) continue;
				if (pad_face_type[f] == 1 && pad_face_type[opp] == 1) {
					pad_face_type[f] = -1;
					pad_face_type[opp] = -1;
				}
			}
			for (auto f : hex.iter_facets()) pad_face[f] = (pad_face_type[f] != -1);
		}


		// DO THE PADDING

		HexPad padder(hex);
		padder.apply(pad_face,false, geom_optim);

		hex.connect();
		halfedge2refine.resize(hex.nfacets() * 4, false);

		for (auto h : hex.iter_halfedges()) {
			if (!pad_face[h.facet()]) continue;
			auto opp = h.opposite_c();
			um_assert(opp.active());
			// edges that are split in this step are not to be refined again
			halfedge2refine[h.opposite_f().next()] = false;
			halfedge2refine[h.opposite_f().prev()] = false;
			halfedge2refine[opp.opposite_f().next()] = false;
			halfedge2refine[opp.opposite_f().prev()] = false;

			// duplicated edges in other directions need to propagate their "to refine" attribute
			halfedge2refine[opp] = halfedge2refine[h];
			halfedge2refine[opp.opposite_f().next().next()] = halfedge2refine[h];
			halfedge2refine[opp.opposite_f().next().next().opposite_f()] = halfedge2refine[h];
		}
		// --- propagate halfedge2refine
		for (auto h : hex.iter_halfedges()) if (halfedge2refine[h] || halfedge2refine[h.opposite_f()])
			for (auto cir : h.iter_CCW_around_edge())
				halfedge2refine[cir] = true;
		//check_halfedge2refine();
	}
}


