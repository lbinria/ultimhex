#pragma once

#include <ultimaille/all.h>
#include <geogram/mesh/mesh_io.h> 

using namespace UM;

namespace algo {

	int tag(Triangle3 &t);

	void naive_tag(Tetrahedra &m, CellFacetAttribute<int> &tag_attr);
	void naive_tag(Triangles &m, FacetAttribute<int> &tag_attr);
	void naive_tag(std::vector<Surface::Facet> facets, FacetAttribute<int> &tag_attr, bool constraints[3]);
	void naive_tag(std::vector<Surface::Facet> facets, GEO::Attribute<GEO::signed_index_t> &flag_attr, bool constraints[3]);
}