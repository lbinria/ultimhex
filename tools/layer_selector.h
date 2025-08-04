#include <vector>
#include "../context.h"

// TODO move here get_facets / get_cells, content of helpers ...
// Use two LayerSelection in LayerSelector (hovered and selected)
struct LayerSelection {

};

struct HalfedgeBag {
	
	HalfedgeBag(Volume &m, std::vector<int> halfedges) : m(m), halfedges(halfedges) {}
	HalfedgeBag(const HalfedgeBag &other) : m(m), halfedges(other.halfedges) {}

	// TODO improve by returning a std::vector<Volume::Facet>
	std::vector<int> get_facets_indexes() {
		std::unordered_set<int> unique_facets;

		for (auto hi : halfedges) {
			Volume::Halfedge h(m, hi);
			unique_facets.insert(h.prev().opposite_f().facet());
		}

		return std::vector(unique_facets.begin(), unique_facets.end());
	}

	std::vector<std::pair<int, int >> get_facets_indexes2() {
		std::set<std::pair<int, int>> unique_facets;

		for (auto hi : halfedges) {
			Volume::Halfedge h(m, hi);
			unique_facets.insert({h.prev().opposite_f().facet(), h.next().opposite_f().facet()});
		}

		return std::vector(unique_facets.begin(), unique_facets.end());
	}

	std::vector<std::array<int, 3>> get_cells_and_facets() {
		std::set<std::array<int, 3>> unique_facets;

		for (auto hi : halfedges) {
			Volume::Halfedge h(m, hi);
			unique_facets.insert({h.cell(), h.prev().opposite_f().facet(), h.next().opposite_f().facet()});
		}

		return std::vector(unique_facets.begin(), unique_facets.end());
	}

	std::vector<Volume::Facet> get_facets(Volume &m) {
		auto facet_indexes = get_facets_indexes();
		std::vector<Volume::Facet> facets;
		for (auto fi : facet_indexes) {
			Volume::Facet f(m, fi);
			facets.push_back(std::move(f));
		}
		return facets;
	}

	std::vector<int> get_cells_indexes() {
		std::unordered_set<int> unique_cells;

		for (auto hi : halfedges) {
			Volume::Halfedge h(m, hi);
			unique_cells.insert(h.cell());
		}

		return std::vector(unique_cells.begin(), unique_cells.end());
	}

	private:
	Volume &m;
	std::vector<int> halfedges;
};

struct FacetBag {

	FacetBag(Volume &m, std::vector<int> facets) : m(m), facets(facets) {
		compute_edges();
	}
	FacetBag(const FacetBag &other) : m(other.m), facets(other.facets) {
		compute_edges();
	}

	std::vector<int> get_edges() const {

		return edges;
		// std::unordered_set<int> unique_edges;
		// std::unordered_set<int> unique_edges_set;

		// for (auto fi : facets) {
		// 	auto f = Volume::Facet(m, fi);
		// 	for (auto h : f.iter_halfedges()) {
		// 		auto oh = h.opposite_f().next();
		// 		if (unique_edges_set.find(oh.from()) == unique_edges_set.end()) {
		// 			unique_edges_set.insert(oh.from());
		// 			unique_edges.insert(oh);
		// 		}
		// 	}
		// }

		// return std::vector(unique_edges.begin(), unique_edges.end());
	}

	// int get_edge(int h) const {
	// 	return halfedge_2_edge[h];
	// }

	private:
	Volume &m;
	std::vector<int> facets;
	std::vector<int> edges;
	// std::vector<int> halfedge_2_edge;

	void compute_edges() {

		edges.clear();

		// std::vector<int> edge_2_halfedge(m.nverts(), -1);
		std::vector<bool> is_process(m.nverts(), false);

		for (auto fi : facets) {
			auto f = Volume::Facet(m, fi);
			for (auto h : f.iter_halfedges()) {
				auto oh = h.opposite_f().next();
				if (!is_process[oh.from()]) {
					// halfedge_2_edge[oh] = edges.size();
					// edge_2_halfedge[edges.size()] = oh;
					is_process[oh.from()] = true;
					edges.push_back(oh);
				} 
				// else {
					// halfedge_2_edge[oh] = is_process[oh.from()];
				// }
			}
		}

	}

	// void compute_edges() {

	// 	halfedge_2_edge.resize(m.ncells() * 24, -1);
	// 	std::unordered_set<int> unique_edges;
	// 	std::unordered_set<int> unique_edges_set;

	// 	for (auto fi : facets) {
	// 		auto f = Volume::Facet(m, fi);
	// 		for (auto h : f.iter_halfedges()) {
	// 			auto oh = h.opposite_f().next();
	// 			if (unique_edges_set.find(oh.from()) == unique_edges_set.end()) {
	// 				unique_edges_set.insert(oh.from());
	// 				unique_edges.insert(oh);
	// 			}
	// 			halfedge_2_edge[h] = oh;
	// 		}
	// 	}

	// 	edges = std::vector(unique_edges.begin(), unique_edges.end());

	// }

};

struct LayerSelector {

	LayerSelector(Context &ctx) : ctx(ctx) {}

	void init();
	void hover();
	void select();
	void highlight_hovered_cells();
	void highlight_selected_cells();
	void clear();

	const std::vector<int> get_selected_halfedges() { return selected_halfedges; }
	const std::vector<int>& get_hovered_halfedges() { return hovered_halfedges; }
	const int get_selected_layer() { return selected_layer; }


	private:
	Context &ctx;
	std::vector<int> layers;

	int hovered_h = -1;

	std::vector<int> hovered_halfedges;
	std::vector<int> selected_halfedges;

	int selected_layer = -1;

};