#include "pointset_in_surface.h"

PointSetEmbedding::PointSetEmbedding(PointSet& pts) : pts(pts), feature(tri_emb), dim(pts, -1), id(pts, -1) {}

void PointSetEmbedding::init_from_triangles(Triangles& m, CornerAttribute<bool>* p_feature) {
    ToolBox(tri_emb).copy_from(m, true);
    tri_emb.connect();
    pts_emb = tri_emb.points;
    pl_emb.points = tri_emb.points;
    for (auto h : tri_emb.iter_halfedges()) feature[h] = (p_feature != NULL && (*p_feature)[h]);
    for (auto h : tri_emb.iter_halfedges()) if (feature[h]) {
        int e = pl_emb.create_edges(1);
        pl_emb.vert(e, 0) = h.from();
        pl_emb.vert(e, 1) = h.to();
    }
    pl_emb.connect();
}

void PointSetEmbedding::set_embedding(int v, int p_dim, int p_id) {
    dim[v] = p_dim;
    id[v] = p_id;
}

void PointSetEmbedding::move_toward_point(int v, vec3 objective) {
    vec3& P = pts[v];
    vec3 dir = objective - P;
    double step = dir.norm();
    dir.normalize();

    if (dim[v] == 0) {
        P = pts_emb[id[v]];
    }
    if (dim[v] == 1) {
        PointOnPolyLine pt(PolyLine::Edge(pl_emb, id[v]), P);
        pt.walk(dir, step);
        P = pt;
        id[v] = pt.e;
    }
    if (dim[v] == 2) {
        PointOnTriangles pt(Surface::Facet(tri_emb, id[v]), P);
        pt.walk(dir, step, &feature);
        P = pt;
        id[v] = pt.f;
    }
}

void PointSetEmbedding::move_toward(PointAttribute<vec3>& objective) {
    FOR(v, pts.size()) move_toward_point(v, objective[v]);
}

void PointSetEmbedding::project_point(int v) {
    if (dim[v] == 0) {
        pts[v] = pts_emb[id[v]];
    }
    if (dim[v] == 1) {
        PointOnPolyLine pt(PolyLine::Edge(pl_emb, id[v]), pts[v]);
        pts[v] = pt.project();
    }
    if (dim[v] == 2) {
        PointOnTriangles pt(Surface::Facet(tri_emb, id[v]), pts[v]);
        pts[v] = pt.project();
    }
}

 void PointSetEmbedding::project() {
    FOR(v, pts.size()) project_point(v);
}

 std::array<vec3, 3> PointSetEmbedding::constrained_direction(int v) {
    if (dim[v] == 0) return { vec3(1,0,0),vec3(0,1,0),vec3(0,0,1) };
    if (dim[v] == 1) {
        auto e = PolyLine::Edge(pl_emb, id[v]);
        vec3 dir = Segment3(e).vector().normalized();



        vec3 x(1, 0, 0);
        if (std::abs(dir * x) > .8) x = vec3(0, 1, 0);
        vec3 y = cross(dir, x).normalized();
        x = cross(dir, y).normalized();
        return { x,y,vec3(0,0,0) };
    }
    if (dim[v] == 2) return { Triangle3(Surface::Facet(tri_emb,id[v])).normal(),vec3(0,0,0),vec3(0,0,0) };
    return { vec3(0,0,0),vec3(0,0,0),vec3(0,0,0) };
}

 void PointSetEmbedding::show_pts() {
    Drop(pts, id).apply();
}

 void PointSetEmbedding::show_emb() {

    DropPointSet(pts_emb).apply("pts_emb");

    EdgeAttribute<int> eid(pl_emb); FOR(e, pl_emb.nedges()) eid[e] = e;
    Drop(pl_emb, eid).apply_wireframe("pl_emb");

    FacetAttribute<int> fid(tri_emb); FOR(f, tri_emb.nfacets()) fid[f] = f;
    Drop(tri_emb, fid).apply("tri_emb");
}

void test_PointSetEmbedding(Triangles& tri) {


    BBox3 box; box.add(-vec3(1, 1, 1)); box.add(vec3(1, 1, 1));
    //ToolBox(tri.points).normalize(box);
    CornerAttribute<int> feature_id(tri);
    FeatureEdgeDetector(tri).dihedral_angle().threshold().remove_small_features().remove_small_features().remove_small_features().apply(feature_id, false);
    Drop(tri, feature_id)._wireframe(true).apply_half_edge("features");

    CornerAttribute<bool> feature(tri);
    for (auto h : tri.iter_halfedges()) feature[h] = (feature_id[h] >= 0);


    PointSetEmbedding emb(tri.points);
    emb.init_from_triangles(tri, &feature);

    PointAttribute<int> nfeats(tri, 0);

    for (auto e : emb.pl_emb.iter_edges()) nfeats[e.from()]++;

    //for (auto h : tri.iter_halfedges()) if (feature[h]) nfeats[h.from()]++; // WARNING boundary is not well managed here 


    for (auto h : tri.iter_halfedges()) {
        Surface::Vertex v = h.from();
        if (emb.dim[v] != -1) continue;
        if (nfeats[v] == 1 || nfeats[v] > 2) { emb.dim[v] = 0; emb.id[v] = v; }
        if (nfeats[v] == 0) { emb.dim[v] = 2; emb.id[v] = h.facet(); }
    }

    for (auto e : emb.pl_emb.iter_edges()) {
        auto v = e.from();
        int val = 0;
        for (auto cir : v.iter_edges()) val++;
        if (val == 2) { emb.dim[v] = 1; emb.id[v] = e; }
    }


    FOR(it, 10) {
        PointAttribute<vec3> smoothed_pos(emb.pts);
        {
            LeastSquares ls(3 * tri.nverts());
            for (auto v : tri.iter_vertices()) {
                auto constraints = emb.constrained_direction(v);//point_handle.constrained_directions();
                for (vec3 c : constraints) {
                    if (c.norm2() < 1e-20) continue;
                    ls.add_to_energy(100. * (c[0] * X(3 * v + 0) + c[1] * X(3 * v + 1) + c[2] * X(3 * v + 2) - c * v.pos()));
                }
            }

            //for (auto h : tri.iter_halfedges())FOR(d, 3)
            //    ls.add_to_energy(X(3 * h.from() + d) - X(3 * h.to() + d));
            for (auto v : tri.iter_vertices()) FOR(d, 3) {
                LinExpr line;
                for (auto h : v.iter_halfedges()) {
                    line += X(3 * h.from() + d) - X(3 * h.to() + d);
                }
                //if (d == 0) line += LinExpr(10);
                ls.add_to_energy(line);
            }


            ls.solve();
            for (auto v : tri.iter_vertices()) FOR(d, 3) smoothed_pos[v][d] = ls.X[3 * v + d];// - v.pos()[d];
        }


        emb.move_toward(smoothed_pos);
        DropSurface(tri).apply("advect");
    }
    DropSurface(tri).apply("tri");



}
