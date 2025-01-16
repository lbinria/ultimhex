#include "helpers.h"
#include <set>
#include <map>
#include <optional>
#include "geom_ultimaille_binding.h"

using namespace UM;

namespace helpers {

	bool left_halfedge(Volume::Halfedge &h, Volume::Halfedge &out_h) {

		auto opp_c = h.next().next().opposite_f().opposite_c();

		if (opp_c.active()) {
			out_h = opp_c.opposite_f();
			return true;
		} else {
			out_h = -1;
			return false;
		}
	}

	bool right_halfedge(Volume::Halfedge &h, Volume::Halfedge &out_h) {
		auto opp_c = h.opposite_f().opposite_c();

		if (opp_c.active()) {
			out_h = opp_c.opposite_f().next().next();
			return true;
		} else {
			out_h = -1;
			return false;
		}
	}

	bool down_halfedge(Volume::Halfedge &h, Volume::Halfedge &out_h) {
		auto opp_c = h.opposite_f().next().next().opposite_f().opposite_c();

		if (opp_c.active()) {
			out_h = opp_c;
			return true;
		} else {
			out_h = -1;
			return false;
		}
	}

	bool up_halfedge(Volume::Halfedge &h, Volume::Halfedge &out_h) {
		auto opp_c = h.opposite_c();

		if (opp_c.active()) {
			out_h = opp_c.opposite_f().next().next().opposite_f();
			return true;
		} else {
			out_h = -1;
			return false;
		}
	}

	bool front_halfedge(Volume::Halfedge &h, Volume::Halfedge &out_h) {
		auto opp_c = h.next().opposite_f().opposite_c();

		if (opp_c.active()) {
			out_h = opp_c.opposite_f().next();
			return true;
		} else {
			out_h = -1;
			return false;
		}
	}

	bool back_halfedge(Volume::Halfedge &h, Volume::Halfedge &out_h) {
		auto opp_c = h.prev().opposite_f().opposite_c();

		if (opp_c.active()) {
			out_h = opp_c.opposite_f().prev();
			return true;
		} else {
			out_h = -1;
			return false;
		}
	}

	int get_facets_layers(UM::Hexahedra &hex, CellFacetAttribute<int> &layer) {

		DisjointSet ds(hex.nfacets());

		for (auto h : hex.iter_halfedges()) {
			Volume::Halfedge l_h(hex, -1);
			Volume::Halfedge b_h(hex, -1);
			auto opp_c = h.next().opposite_f().opposite_c();
			if (opp_c.active()) {
				ds.merge(h.facet(), opp_c.opposite_f().next().facet());
			}
			opp_c = h.opposite_f().opposite_c();
			if (opp_c.active()) {
				ds.merge(h.facet(), opp_c.opposite_f().next().next().facet());
			}
		}

		return ds.get_sets_id(layer.ptr->data);
	}

	std::vector<int> get_cells_layer(Hexahedra &hex) {
		
		DisjointSet ds(hex.ncells());

		for (auto h : hex.iter_halfedges()) {
			auto opp = h.next().next().opposite_f().opposite_c();
			if (opp.active()) {
				ds.merge(h, opp.opposite_f());
			}
			opp = h.opposite_f().next().next().opposite_f().opposite_c();
			if (opp.active()) {
				ds.merge(h, opp);
			}
		}

		std::vector<int> layers;
		ds.get_sets_id(layers);
		return layers;
	}

	std::vector<int> get_h_layers(Hexahedra &hex) {
		
		DisjointSet ds(24 * hex.ncells());

		for (auto h : hex.iter_halfedges()) {

			auto opp = h.next().next().opposite_f().opposite_c();

			if (opp.active())
				ds.merge(h, opp.opposite_f());
			
			opp = h.opposite_c();

			if (opp.active())
				ds.merge(h, opp.opposite_f().next().next().opposite_f());
		}

		std::vector<int> layers;
		ds.get_sets_id(layers);
		return layers;
	}

	int get_layers(Hexahedra &hex, EdgeGraph &eg, EdgeAttribute<int> &layer) {
		
		DisjointSet ds(eg.nedges());

		for (auto h : hex.iter_halfedges()) {
			ds.merge(eg.edge_from_halfedge(h), eg.edge_from_halfedge(h.next().next().opposite_f()));
		}

		return ds.get_sets_id(layer.ptr->data);
	}

	void get_layer_graph(Hexahedra &hex, EdgeGraph &eg, EdgeAttribute<int> &layer, int nlayers, EdgeAttribute<int> &next_eg, EdgeAttribute<int> &prev_eg, PolyLine &layer_graph) {
		
		for (auto h : hex.iter_halfedges()) {

			auto e0 = eg.edge_from_halfedge(h);
			auto next = h.next().opposite_f().opposite_c();

			if (!next.active()) {
				next_eg[e0] = -1;
				continue;
			}

			next = next.opposite_f().next();
			auto e1 = eg.edge_from_halfedge(next);

			if (next_eg[e0] == -2) 
				next_eg[e0] = e1;
			// Many next => singus
			else if (next_eg[e0] != e1) 
				next_eg[e0] = -1;
		}

		// Set previous from next
		for (auto e : eg.iter_edges()) 
			if (next_eg[e] >= 0) 
				prev_eg[next_eg[e]] = e;

		// ---------------------------

		// Create polyline
		{
			layer_graph.points.create_points(nlayers);
			layer_graph.connect();
			PointAttribute<int> size(layer_graph, 0);
			
			for (auto e : eg.iter_edges()) { 
				size[layer[e]]++; 
				// layer_graph.points[layer[e]] += .8 * e.from().pos() + .2 * e.to().pos();
				layer_graph.points[layer[e]] += .5 * e.from().pos() + .5 * e.to().pos();
			}

			for (auto v : layer_graph.iter_vertices()) 
				v.pos() /= double(size[v]);

			
			for (auto e : eg.iter_edges()) if (next_eg[e] >= 0) {
				PolyLine::Vertex v0(layer_graph, layer[e]);
				PolyLine::Vertex v1(layer_graph, layer[next_eg[e]]);
			
				bool already_in = false;
				for (auto e : v0.iter_edges())
					already_in = already_in || (e.to() == v1);
			
				if (!already_in) 
					layer_graph.conn->create_edge(v0, v1);
			}
		}


		std::vector<bool> to_kill(layer_graph.nedges(), false);

		for (auto e : eg.iter_edges()) {
			PolyLine::Vertex v0(layer_graph, layer[e]);

			if (next_eg[e] < 0)
				for (auto e_l : v0.iter_edges())
					to_kill[e_l] = true;

			if (prev_eg[e] < 0)
				for (auto e_l : layer_graph.iter_edges())
					if (e_l.to() == layer[e])
						to_kill[e_l] = true;
		}

		layer_graph.delete_edges(to_kill);
		layer_graph.connect();
	}

	std::vector<UM::vec3> get_layer_stack_visu(Hexahedra &hex, Volume::Halfedge &selected_he) {
		assert(hex.connected());

		// Get layers
		EdgeGraph eg(hex);
		EdgeAttribute<int> layer(eg);
		int nlayers = get_layers(hex, eg, layer);

		// Get layer graph
		PolyLine layer_graph;
		EdgeAttribute<int> next_eg(eg, -2);
		EdgeAttribute<int> prev_eg(eg, -1);
		get_layer_graph(hex, eg, layer, nlayers, next_eg, prev_eg, layer_graph);

		// Group stack of layers together
		PointAttribute<int> layer_stack(layer_graph);
		{
			DisjointSet ds(layer_graph.nverts());
			for (auto e : layer_graph.iter_edges())
				ds.merge(e.from(), e.to());

			ds.get_sets_id(layer_stack.ptr->data);
		}


		auto selected_edge = eg.edge_from_halfedge(selected_he);
		int stack = layer_stack[layer[selected_edge]];


		std::vector<UM::vec3> points;
		for (auto seed : layer_graph.iter_edges()) {
			if (layer_stack[seed.from()] == stack) {
				points.push_back(seed.from().pos());
			}
		}

		return points;
	}

	// Get layer stack cells
	std::vector<std::pair<int, int>> get_layer_stack(Hexahedra &hex, Volume::Halfedge &selected_he) {

		assert(hex.connected());

		// Get layers
		EdgeGraph eg(hex);
		EdgeAttribute<int> layer(eg);
		int nlayers = get_layers(hex, eg, layer);


		// Get layer graph
		PolyLine layer_graph;
		EdgeAttribute<int> next_eg(eg, -2);
		EdgeAttribute<int> prev_eg(eg, -1);
		get_layer_graph(hex, eg, layer, nlayers, next_eg, prev_eg, layer_graph);

		// Group stack of layers together
		PointAttribute<int> layer_stack(layer_graph);
		{
			DisjointSet ds(layer_graph.nverts());
			for (auto e : layer_graph.iter_edges())
				ds.merge(e.from(), e.to());

			ds.get_sets_id(layer_stack.ptr->data);
		}

		auto selected_edge = eg.edge_from_halfedge(selected_he);
		int stack = layer_stack[layer[selected_edge]];
		
		CellFacetAttribute<bool> fa(hex, false);

		// // Extract edges from stack
		std::vector<std::pair<int,int>> cell_facets;
		// for (auto e : eg.iter_edges()) {
		// 	auto h = eg.halfedge_from_edge(e);
		// 	auto c = h.cell();
		// 	auto f = h.prev().opposite_f().facet();
		// 	fa[f] = layer_stack[layer[e]];
		// 	ea[e] = layer_stack[layer[e]];

		// 	if (layer_stack[layer[e]] == stack) {
				
		// 		auto gf = um_bindings::geo_facet_index_from_um_facet_index(f, 6);
		// 		auto lf = um_bindings::geo_local_cell_facet_index_from_facet(gf);
		// 		cell_facets.push_back({c, lf});
		// 	}
		// }

		for (auto f : hex.iter_facets()) {

			auto e = eg.edge_from_halfedge(f.halfedge(0).opposite_f().next());

			if (layer_stack[layer[e]] == stack) {
				fa[f] = true;
				// fa[f] = layer_stack[layer[e]];
				auto c = f.cell();
				auto gf = um_bindings::geo_facet_index_from_um_facet_index(f, 6);
				auto lf = um_bindings::geo_local_cell_facet_index_from_facet(gf);
				cell_facets.push_back({c, lf});
			}
		}


		UM::Quads q_out;
		for (auto f : hex.iter_facets()) {
			// if (fa[f] < 0)
			if (!fa[f])
				continue;

			int v_off = q_out.points.create_points(4);
			q_out.points[v_off] = f.vertex(0).pos();
			q_out.points[v_off + 1] = f.vertex(1).pos();
			q_out.points[v_off + 2] = f.vertex(2).pos();
			q_out.points[v_off + 3] = f.vertex(3).pos();

			int f_off = q_out.create_facets(1);
			q_out.vert(f_off, 0) = v_off;
			q_out.vert(f_off, 1) = v_off + 1;
			q_out.vert(f_off, 2) = v_off + 2;
			q_out.vert(f_off, 3) = v_off + 3;

		};
		write_by_extension("q_out.geogram", q_out, {{}, {/*{"f", fa.ptr}*/}, {}});





		return cell_facets;
	}

	void redefine_stack_layers(Hexahedra &hex, Volume::Halfedge &selected_he, int final_height) {
		assert(hex.connected());

		std::vector<bool> cells_to_kill(hex.ncells());
		std::vector<vec3> npts;
		std::vector<int> new_hexes;

		{
			// Get layers
			EdgeGraph eg(hex);
			EdgeAttribute<int> layer(eg);
			int nlayers = get_layers(hex, eg, layer);

			// Get layer graph
			PolyLine layer_graph;
			EdgeAttribute<int> next_eg(eg, -2);
			EdgeAttribute<int> prev_eg(eg, -1);
			get_layer_graph(hex, eg, layer, nlayers, next_eg, prev_eg, layer_graph);

			// Group stack of layers together
			PointAttribute<int> layer_stack(layer_graph);
			{
				DisjointSet ds(layer_graph.nverts());
				for (auto e : layer_graph.iter_edges())
					ds.merge(e.from(), e.to());

				ds.get_sets_id(layer_stack.ptr->data);
			}

			// Get selected edge from selected halfedge
			// Retrieve selected stack layer
			auto selected_edge = eg.edge_from_halfedge(selected_he);
			int stack = layer_stack[layer[selected_edge]];


			CellFacetAttribute<int> start(hex, -1);
			CellFacetAttribute<bool> fa(hex, false);

			int n_start = 0;


			for (auto f : hex.iter_facets()) {

				// auto e = eg.edge_from_halfedge(f.halfedge(0).opposite_f().next());

				for (auto h : f.iter_halfedges()) {
					auto e = eg.edge_from_halfedge(h.opposite_f().next());

					if (layer_stack[layer[e]] == stack) {
						fa[f] = true;
					}

					// Should fix here !
					if (layer_stack[layer[e]] == stack && (prev_eg[e] < 0 || layer_stack[layer[prev_eg[e]]] != stack)) {
						start[f] = 1;
						break;
					}
				}


			}


			for (auto f : hex.iter_facets()) {
				if (start[f] > 0)
					n_start++;
			}

			std::cout << "n start: " << n_start << std::endl;

			UM::Quads q_out;
			FacetAttribute<int> start2(q_out, 0);
			for (auto f : hex.iter_facets()) {
				// if (fa[f] < 0)
				if (!fa[f])
					continue;

				
				int v_off = q_out.points.create_points(4);
				q_out.points[v_off] = f.vertex(0).pos();
				q_out.points[v_off + 1] = f.vertex(1).pos();
				q_out.points[v_off + 2] = f.vertex(2).pos();
				q_out.points[v_off + 3] = f.vertex(3).pos();

				int f_off = q_out.create_facets(1);
				q_out.vert(f_off, 0) = v_off;
				q_out.vert(f_off, 1) = v_off + 1;
				q_out.vert(f_off, 2) = v_off + 2;
				q_out.vert(f_off, 3) = v_off + 3;

				start2[f_off] = start[f];

			};
			write_by_extension("q_out.geogram", q_out, {{}, {{"f", start2.ptr}}, {}});

			// Map origin vertex id of poly to it's next vertices 
			std::map<int, std::vector<int> > extrusion;

			{
				for (auto f : hex.iter_facets()) {

					// Search for a starting facet
					if (start[f] <= 0)
						continue;

					for (auto h : f.iter_halfedges()) {

						auto e = eg.edge_from_halfedge(h.opposite_f().next());

						auto org = e.from();
						
						if (extrusion.find(org) != extrusion.end())
							continue;

						extrusion[org] = std::vector<int>(1, org);

						// Search for other extremity
						std::vector<double> segment_lengths;
						std::vector<UM::vec3> segment_points;
						
						segment_points.push_back(org.pos());
						segment_lengths.push_back((e.to().pos() - e.from().pos()).norm());

						while (next_eg[e] > 0 && layer_stack[layer[next_eg[e]]] == stack) {
							auto next_e = PolyLine::Edge(eg, next_eg[e]);
							e = next_eg[e];
							segment_points.push_back(e.from().pos());
							double l = (e.to().pos() - e.from().pos()).norm();
							segment_lengths.push_back(l);
						}
						auto dest = e.to();
						
						segment_points.push_back(e.to().pos());

						// Compute size of polyline
						double poly_len = 0;
						for (auto l : segment_lengths)
							poly_len += l;

						double segment_len_target = poly_len / final_height;


						int offv =  npts.size();
						npts.resize(offv + final_height - 1);
						for (int lv = 0; lv < final_height - 1; lv++) {
							// Push new vertex id to polyline (combinatoric)
							extrusion[org].push_back(hex.nverts() + offv + lv);

							// Compute geometric pos
							double t = segment_len_target * (lv + 1);
							// Search on which source segment this point is located !
							double cur_t = t;
							double i = 0;
							for (auto length : segment_lengths) {
								cur_t -= length;
								if (cur_t < 0) {
									cur_t += length;
									break;
								}

								++i;
							}

							// Make interpolation between node of choosen segment
							auto dir = segment_points[i + 1] - segment_points[i];
							auto pos = segment_points[i] + dir * (cur_t / segment_lengths[i]);

							// Set pos (geometry)
							npts[offv + lv] = pos;

						}

						std::cout << "npts resize: " << npts.size() << std::endl;

						extrusion[org].push_back(dest);

					}
				}

			}

			std::cout << "extrusion.size() = " << extrusion.size() << std::endl;


			// {
			// 	for (auto seed : eg.iter_edges()) {

			// 		// auto seed_h = eg.halfedge_from_edge(seed);
			// 		// auto f = seed_h.prev().opposite_f();
			// 		// if (start[f] <= 0)
			// 		// 	continue;

			// 		if (!(layer_stack[layer[seed]] == stack && (prev_eg[seed] < 0 || layer_stack[layer[prev_eg[seed]]] != stack)))
			// 			continue;

			// 		// // Search for edge seed + Search for seed of requested stack
			// 		// if (prev_eg[seed] < 0 && layer_stack[layer[seed]] == stack) {

			// 			n_origin++;

			// 			auto e = seed;

			// 			auto org = e.from();
			// 			extrusion[org] = std::vector<int>(1, org);

			// 			// Search for other extremity
			// 			while (next_eg[e] > 0 && layer_stack[layer[next_eg[e]]] == stack) {
			// 				e = next_eg[e];
			// 			}
			// 			auto dest = e.to();
			// 			std::cout << std::endl;

			// 			int offv =  npts.size();
			// 			npts.resize(offv + final_height - 1);
			// 			for (int lv = 0; lv < final_height - 1; lv++) {
			// 				// Push new vertex id to polyline (combinatoric)
			// 				extrusion[org].push_back(hex.nverts() + offv + lv);
			// 				// Set pos (geometry)

			// 				double c = double(lv + 1) / double(final_height);
			// 				npts[offv + lv] = (1. - c) * org.pos() + c * dest.pos();

			// 			}
			// 			extrusion[org].push_back(dest);
			// 		}

			// }

			// std::cout << "extrusion size: " << extrusion.size() << std::endl;

			Quads qq;
			qq.points.create_points(8);

			for (auto f : hex.iter_facets()) {
				if (start[f] <= 0)
					continue;

				// Assert that each halfedge has an extrusion
				for (auto h : f.iter_halfedges()) 
					um_assert(extrusion.find(h.from()) != extrusion.end());
				
				// std::vector<int>& stack0 = extrusion[f.halfedge(0).from()];
				// std::vector<int>& stack1 = extrusion[f.halfedge(1).from()];
				// std::vector<int>& stack2 = extrusion[f.halfedge(3).from()];
				// std::vector<int>& stack3 = extrusion[f.halfedge(2).from()];

				std::vector<int>& stack0 = extrusion[f.halfedge(1).from()];
				std::vector<int>& stack1 = extrusion[f.halfedge(0).from()];
				std::vector<int>& stack2 = extrusion[f.halfedge(2).from()];
				std::vector<int>& stack3 = extrusion[f.halfedge(3).from()];
				
				for (int l = 0; l < final_height; l++) { 
					for (int dl = 0; dl < 2; dl++){
						new_hexes.push_back(stack0[l+dl]);
						new_hexes.push_back(stack1[l+dl]);
						new_hexes.push_back(stack2[l+dl]);
						new_hexes.push_back(stack3[l+dl]);
					}
				}

				// auto p0 =  extrusion[f.halfedge(0).from()][0];
				// auto p1 =  extrusion[f.halfedge(1).from()][0];
				// auto p3 =  extrusion[f.halfedge(3).from()][0];
				// auto p2 =  extrusion[f.halfedge(2).from()][0];
				// auto p4 =  extrusion[f.halfedge(0).from()][1];
				// auto p5 =  extrusion[f.halfedge(1).from()][1];
				// auto p7 =  extrusion[f.halfedge(3).from()][1];
				// auto p6 =  extrusion[f.halfedge(2).from()][1];

				// qq.points[0] = npts[p0];
				// qq.points[1] = npts[p1];
				// qq.points[2] = npts[p2];
				// qq.points[3] = npts[p3];
				// qq.points[4] = npts[p4];
				// qq.points[5] = npts[p5];
				// qq.points[6] = npts[p6];
				// qq.points[7] = npts[p7];
				
				// qq.create_facets(3);
				// qq.vert(0, 0) = p0;
				// qq.vert(0, 1) = p1;
				// qq.vert(0, 2) = p2;
				// qq.vert(0, 3) = p3;

				// qq.vert(1, 0) = p4;
				// qq.vert(1, 1) = p5;
				// qq.vert(1, 2) = p6;
				// qq.vert(1, 3) = p7;
				// break;
			}

			
			for (auto h : hex.iter_halfedges()) 
				if (layer_stack[layer[eg.edge_from_halfedge(h)]] == stack) 
					cells_to_kill[h.cell()] = true;
		}

		hex.delete_cells(cells_to_kill);

		std::cout << "create npts: " << npts.size() << std::endl;

		int offv = hex.points.create_points(npts.size());
		for (int i = 0; i < npts.size(); i++) 
			hex.points[offv + i] = npts[i];

		int offc = hex.create_cells(new_hexes.size() / 8);

		for (int c = 0; c < new_hexes.size(); c++) 
			hex.vert(offc + (c / 8), c % 8) = new_hexes[c];
		
		hex.delete_isolated_vertices();
	}

	// Works !
	void layer_from_halfedge(UM::Hexahedra &hex, int he, std::function<void(UM::Volume::Halfedge&)> f) {
		// Must be connected
		assert(hex.connected());

		std::vector<bool> visited(24 * hex.ncells(), false);
		std::queue<int> q;
		q.push(he);
		visited[he] = true;

		while (!q.empty()) {

			auto he_idx = q.front();
			auto he = Volume::Halfedge(hex, he_idx);
			q.pop();

			// Front unfold
			auto opp_c = he.next().next().opposite_f().opposite_c();
			if (opp_c.active()) {
				auto nxt_he = opp_c.opposite_f();
				if (!visited[nxt_he]) {
					q.push(nxt_he);
					visited[nxt_he] = true;
				}
			}
			// Back unfold
			opp_c = he.opposite_f().opposite_c();
			if (opp_c.active()) {
				auto nxt_he = opp_c.opposite_f().next().next();
				if (!visited[nxt_he]) {
					q.push(nxt_he);
					visited[nxt_he] = true;
				}
			}
			// Down unfold
			opp_c = he.opposite_f().next().next().opposite_f().opposite_c();
			if (opp_c.active()) {
				auto nxt_he = opp_c;
				if (!visited[nxt_he]) {
					q.push(nxt_he);
					visited[nxt_he] = true;
				}
			}
			// Up unfold
			opp_c = he.opposite_c();
			if (opp_c.active()) {
				auto nxt_he = opp_c.opposite_f().next().next().opposite_f();
				if (!visited[nxt_he]) {
					q.push(nxt_he);
					visited[nxt_he] = true;
				}
			}

			// Call callback function passing current processed halfedge
			f(he);
		}
	}
	




	void collapse(UM::Hexahedra &hex, std::vector<int> layer) {
		// Use edge graph
		// Compute indirection map of vertice to vertice, set vertice to vertice at the begging, but set from() to to() for collapsing

		std::vector<int> indirection_map(hex.nverts());
		std::vector<bool> cells_to_kill(hex.ncells(), false);
		{
			// EdgeGraph eg(hex);
			for (auto v : hex.iter_vertices())
				indirection_map[v] = v;

			for (auto h_idx : layer) {
				Volume::Halfedge h(hex, h_idx);
				// auto e = eg.edge_from_halfedge(h);
				indirection_map[h.from()] = h.to();
				cells_to_kill[h.cell()] = true;
			}
		}

		for (int c = 0; c < hex.ncells(); c++) {
			for (int lv = 0; lv < 8; lv++) {
				hex.vert(c, lv) = indirection_map[hex.vert(c, lv)];
			}
		}

		// Remove degenerate hex
		// hex.disconnect();
		// hex.delete_cells(cells_to_kill);

		// Remove isolated vertices
		hex.delete_isolated_vertices();
		// hex.connect();
	}

	void layer_along(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Facet&)> f) {
		assert(hex.connected());

		std::vector<bool> visited(24 * hex.ncells(), false);
		std::queue<int> q;
		q.push(start_he);
		visited[start_he] = true;


		while (!q.empty()) {

			auto he_idx = q.front();
			auto he = Volume::Halfedge(hex, he_idx);

			q.pop();

			auto opp_c = he.next().opposite_f().opposite_c();

			if (opp_c.active()) {
				auto n_he = opp_c.opposite_f().next();

				if (!visited[n_he]) {
					q.push(n_he);
					visited[n_he] = true;
				}
			} else {
				auto n_he = he.next().opposite_f().next();
				if (!visited[n_he]) {
					q.push(n_he);
					visited[n_he] = true;
				}
			}

			opp_c = he.prev().opposite_f().opposite_c();

			if (opp_c.active()) {
				auto n_he = opp_c.opposite_f().prev();

				if (!visited[n_he]) {
					q.push(n_he);
					visited[n_he] = true;
				}
			} else {
				auto n_he = he.prev().opposite_f().prev();
				if (!visited[n_he]) {
					q.push(n_he);
					visited[n_he] = true;
				}
			}

			// Propagate to front
			opp_c = he.opposite_f().next().next().opposite_f().opposite_c();

			if (opp_c.active()) {
				auto n_he = opp_c;

				if (!visited[n_he]) {
					q.push(n_he);
					visited[n_he] = true;
				}
			} 

			// Process current
			auto facet = he.opposite_f().facet();
			f(facet);
			
		}
	}

	void loop_along(UM::Hexahedra &hex, UM::Volume::Halfedge &start_he, std::function<void(UM::Volume::Halfedge& /* cur_he */, bool/* on_border */)> f) {
		assert(hex.connected());
		std::vector<bool> visited(24 * hex.ncells(), false);

		auto cur_he = start_he;


		while (true) {

			f(cur_he, true);

			auto opp_c = cur_he.next().opposite_f().opposite_c();
			if (!opp_c.active() || cur_he.next().opposite_f().next().facet().on_boundary()) {
				cur_he = cur_he.next().opposite_f().next();

			} else if (opp_c.active()) {
				cur_he = opp_c.opposite_f().next();

				// on border ?
			}

			if (visited[cur_he])
				break;

			visited[cur_he] = true;


		}

	}

}