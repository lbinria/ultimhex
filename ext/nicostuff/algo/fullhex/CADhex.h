#ifndef CADHEX__H__
#define CADHEX__H__
#include <ultimaille/all.h>
#include <framework/trace.h>
#include <fullhex/gp_basic.h>



#include<algo/surface/pointset_in_surface.h>

struct HexCAD{
	HexCAD(Hexahedra& hex) ;
	HexCAD& init_hex_from_uvw(TetBoundary &bound, CornerAttribute<bool>& feature, CellCornerAttribute<vec3>& U) ;
	HexCAD& drop_SJ(std::string name ="");
	HexCAD& verbose_level(int i) { verbose = i; return *this; }


	PointSetEmbedding emb;
	int verbose;
	Hexahedra& hex;
};





#endif