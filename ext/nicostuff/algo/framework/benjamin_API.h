#ifndef BENJAMIN__API__H__
#define BENJAMIN__API__H__

#include <ultimaille/all.h>

using namespace UM;
namespace BenjaminAPI {

	void polycubify(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, Hexahedra& hex, int nhex_wanted, CellFacetAttribute<int> &emb_out);
	void polycubify(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, Hexahedra& hex, int nhex_wanted);
	bool pad(Hexahedra& hex,CellFacetAttribute<bool>& pad_face, CellFacetAttribute<int> &emb_attr);
	void smooth(Hexahedra &hex, CellFacetAttribute<int>&emb_attr, Triangles &tri, FacetAttribute<int> &tri_chart);
	void embeditinit(Triangles &tri, FacetAttribute<int> &tri_chart, Hexahedra &hex, CellFacetAttribute<int> &emb_attr, FacetAttribute<int> &quad_chart, bool gmsh_chart);
	void embeditapply(Hexahedra &hex, CellFacetAttribute<int> &emb_attr, Quads &quad, FacetAttribute<int> &quad_chart, Triangles &tri, FacetAttribute<int> &tri_chart);
	
};

#endif