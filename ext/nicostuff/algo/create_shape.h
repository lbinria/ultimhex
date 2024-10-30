#ifndef CREATE_SHAPE_H___
#define CREATE_SHAPE_H___

#include <ultimaille/all.h>
using namespace UM;

namespace CreatePointSet {
    void from_vector(PointSet& p, const std::vector<vec3>& data);
    void grid1d(PointSet& p, int n);
    void grid2d(PointSet& p, int n);
    void grid3d(PointSet& p, int n);
};
namespace CreatePolyline {
    void grid1d(PolyLine& m, int n);
    void grid2d(PolyLine& m, int n);
    void circle(PolyLine& m, int n);
};
namespace CreateSurface {
    void quad_grid(Quads& m,int n);
    void triangles_grid(Triangles& m,int n);
    void tet(Triangles& m);
    void non_manifold_vertex1(Triangles& m);
    void non_manifold_vertex2(Triangles& m);
    void non_manifold_edge(Triangles& m);
    void triangles_fan_negative_curvature(Triangles& m);
    void polygonal_hut(Polygons& m,double radius = .5, double height = .5, int n = 20);
};
namespace CreateVolume {
    void edge_one_ring(Tetrahedra& m,int n = 6);
    void one_hex(Hexahedra& m);
    void hex_grid(Hexahedra& m,int n);
    void two_wedges(Wedges& m);
    void pyramids_cube(Pyramids& m);
    void tri_grid(Hexahedra& grid, int nb_subdivision);
};

void test_create_shapes();

#endif