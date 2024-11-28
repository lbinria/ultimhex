
#include <ultimaille/all.h>

namespace helpers {

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