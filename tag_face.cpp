#include <ultimaille/all.h>

using namespace UM;

namespace algo {

	int tag(Triangle3 &t) {
		vec3 n = t.normal();
		for (int d = 0; d < 3; d++) {

			if (std::abs(n[d]) < std::abs(n[(d + 1) % 3]) || 
				std::abs(n[d]) < std::abs(n[(d + 2) % 3]))
				continue;

			if (n[d] < 0) 
				return d;
			else 
				return d+3;
		}

		return -1;
	}

	void naive_tag(Tetrahedra &m, CellFacetAttribute<int> &tag_attr) {
		if (!m.connected()) 
			m.connect();

		for (auto f : m.iter_facets()) {
			if (!f.on_boundary())
				continue;
			
			Triangle3 t = f;
			tag_attr[f] = tag(t);
		}
	}

	void naive_tag(Triangles &m, FacetAttribute<int> &tag_attr) {
		for (auto f : m.iter_facets()) {
			Triangle3 t = f;
			tag_attr[f] = tag(t);
		}
	}

}