#ifndef HEX_EDIT__H__
#define	HEX_EDIT__H__

#include <ultimaille/all.h>
#include <framework/trace.h>


#include "toolbox.h"
#include "drop_attribute.h"
#include <algo/volume/hex_smooth.h>



struct HexPad {
	HexPad(Hexahedra& hex);
	void apply(CellFacetAttribute<bool>& pad_face, bool verbose = false, bool geom_optim = true, bool dilate  =false,PointAttribute<vec3>* other_pos = NULL);
	Hexahedra& hex;

};
void refine_marked_cells(Hexahedra& hex, CellAttribute<bool>& need_refine, bool geom_optim = true);


#endif