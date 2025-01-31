#ifndef BENJAMIN__API__H__
#define BENJAMIN__API__H__

#include <ultimaille/all.h>
using namespace UM;
namespace BenjaminAPI {
	void polycubify(Tetrahedra& tet, CellFacetAttribute<int>& tet_flag, Hexahedra& hex, int nhex_wanted);
	bool pad(Hexahedra& hex,CellFacetAttribute<bool>& pad_face);
};

#endif