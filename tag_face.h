#pragma once

#include <ultimaille/all.h>

using namespace UM;

namespace algo {

	int tag(Triangle3 &t);

	void naive_tag(Tetrahedra &m, CellFacetAttribute<int> &tag_attr);
	void naive_tag(Triangles &m, FacetAttribute<int> &tag_attr);

}