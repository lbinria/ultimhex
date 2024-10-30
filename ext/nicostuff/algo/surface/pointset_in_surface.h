#ifndef POIN_SET_IN_SURFACE__H__
#define POIN_SET_IN_SURFACE__H__

#include <surface/feature_curve_detector.h>
#include <basic.h>
#include <drop_mesh.h>
#include <mesh_geom.h>
#include <ultimaille/all.h>


struct PointSetEmbedding {
    PointSetEmbedding(PointSet& pts);

    void init_from_triangles(Triangles& m,CornerAttribute<bool>* p_feature=NULL);

    void set_embedding(int v,int p_dim,int p_id);

    void move_toward_point(int v, vec3 objective);

    void move_toward(PointAttribute<vec3>& objective);


    void project_point(int v);

    void project();


    std::array<vec3, 3> constrained_direction(int v);

    void show_pts ();
    void show_emb ();

    PointSet& pts;
    PointAttribute<int> dim;
    PointAttribute<int> id;

    PointSet pts_emb;
    PolyLine pl_emb;
    Triangles tri_emb;
    CornerAttribute<bool> feature;
};


using namespace UM::Linear;
void test_PointSetEmbedding(Triangles& tri);;

#endif