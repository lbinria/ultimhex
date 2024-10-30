#include <fullhex/hextract.h>
#include <fullhex/gp_basic.h>
#include <volume/hex_select.h>


DualContour::DualContour(Hexahedra& hex) :hex(hex), cond(hex), modified(hex, false) { verbose_level(1); }

DualContour& DualContour::drop_SJ(std::string name) {
    if (verbose>0) ToolBox(hex).drop_scaled_jacobien(name);
    return *this;
}


const vec3 hex_face_to_normal[6] = {
    vec3(-1, 0, 0), vec3(1, 0, 0),
    vec3(0, -1, 0),vec3(0, 1, 0),
    vec3(0, 0, -1),vec3(0, 0, 1)
};

DualContour& DualContour::untangle(double fit_scale ) {
    UntanglerHexTan untangle(hex);
    untangle.apply(*this,fit_scale);
    return* this;
}


DualContour& DualContour::init_hex_from_uvw(Tetrahedra& tet, CellCornerAttribute<vec3>& uvw) {

    double hex_size = ToolBox(tet).ave_edge_size();
    Drop(tet, uvw).apply("uvw");

    // initialize all cubes
    CellAttribute<int> tet_id(hex);
    CellAttribute<vec3> pos(hex);
    DynamicKNN<3> knn(pos.ptr->data);
    CellAttribute<vec3> pos_uvw(hex);
    Trace::step("init all cubes");
    for (auto tetra : tet.iter_cells()) {
        //if (tetra>0) return;
        BBox3 box;
        Tetrahedron tet_uvw(uvw[tetra.corner(0)], uvw[tetra.corner(1)], uvw[tetra.corner(2)], uvw[tetra.corner(3)]);

        if (tet_uvw.volume() < 1e-10) continue;
        FOR(lv, 4) box.add(tet_uvw[lv]);
        box.dilate(1e-20);
        for (int i = int(std::floor(box.min[0])); i < int(std::ceil(box.max[0])); i++)
            for (int j = int(std::floor(box.min[1])); j < int(std::ceil(box.max[1])); j++)
                for (int k = int(std::floor(box.min[2])); k < int(std::ceil(box.max[2])); k++) {
                    // test is_in in texture space
                    vec3 guess_uvw = vec3(i, j, k) + .5 * vec3(1, 1, 1);
                    vec4 bc = tet_uvw.bary_coords(guess_uvw);
                    bool is_in = true;
                    FOR(lv, 4) is_in = is_in && (bc[lv] > -1e-20);
                    if (!is_in) continue;

                    // create the hex centered on guess_uvw
                    vec3 bary_pos(0, 0, 0);
                    FOR(lv, 4) bary_pos += bc[lv] * tetra.vertex(lv).pos();
                    bool already_in = false;


                    auto nn = knn.query(bary_pos);
                    if (nn.size() > 0) {
                        int check_c = nn[0];
                        //FOR(check_c, hex.ncells()) 
                        already_in = already_in || (pos[check_c] - bary_pos).norm() < hex_size / 100.;
                    }
                    if (already_in) continue;

                    int offset_v = hex.points.create_points(8);
                    int offset_c = hex.create_cells(1);
                    pos_uvw[offset_c] = guess_uvw;
                    pos[offset_c] = bary_pos;
                    tet_id[offset_c] = tetra;
                    FOR(lv, 8) hex.vert(offset_c, lv) = offset_v + lv;

                    mat3x3 J_inv = uvw_to_jacobian(tetra, uvw).invert();
                    FOR(di, 2)FOR(dj, 2)FOR(dk, 2) hex.points[offset_v + di + 2 * dj + 4 * dk] = pos[offset_c] + .5 * J_inv * (vec3(di, dj, dk) - .5 * vec3(1, 1, 1));
                }
    }

    if (verbose > 1)Drop(hex, tet_id).apply("tet_id");

    //return *this;

    PolyLine pl;
    EdgeAttribute<int> plrot(pl);
    DisjointSet ds(hex.nverts());



    hex.connect();
    Trace::step("Find links between cubes");
    // connect cube faces to "opposite OR surface point"
    FOR(cur_hex, hex.ncells()) FOR(lf, 6) {


        vec3 dir = hex_face_to_normal[lf];
        vec3 objective_uvw = pos_uvw[cur_hex] + dir;
        Volume::Cell cur_tet(tet, tet_id[cur_hex]);

        int objective_hex = -1;
        vec3 I;
        vec3 G = pos_uvw[cur_hex];
        GPTransitionFunction tf;
        while (true) {
            Tetrahedron tet_uvw(
                tf.apply(uvw[cur_tet.corner(0)]),
                tf.apply(uvw[cur_tet.corner(1)]),
                tf.apply(uvw[cur_tet.corner(2)]),
                tf.apply(uvw[cur_tet.corner(3)])
            );
            if (tet_uvw.volume() < 0) {
                Trace::alert("Cannot cross an tet with negative volume");
                break;
            }

            vec4 bc = tet_uvw.bary_coords(objective_uvw);

            {// check if the destination is reached
                bool is_in = true;
                FOR(lv, 4) is_in = is_in && (bc[lv] > -1e-10);
                if (is_in) {
                    vec3 objective_xyz = vec3(0, 0, 0);
                    FOR(lv, 4) objective_xyz += bc[lv] * cur_tet.vertex(lv).pos();
                    int inner_cur_hex = knn.query(objective_xyz)[0];
                    if ((pos[inner_cur_hex] - objective_xyz).norm() < 1e-2 * hex_size && tet_id[inner_cur_hex] == cur_tet)
                        objective_hex = inner_cur_hex;
                    else std::cerr << "Found a tet that contains the objective, but does not match an hex :( \n\n";
                    break;
                }
            }


            {// objective_uvw is not in cut_tet, need to visit the next one.
                Volume::Facet out_facet(tet, -1);

                for (auto f : cur_tet.iter_facets()) {
                    Triangle3 tr_uvw({
                        tf.apply(uvw[f.corner(0)]),
                        tf.apply(uvw[f.corner(1)]),
                        tf.apply(uvw[f.corner(2)])
                        });
                    vec3 bc_tri;

                    if (!Intersect::triangle_line(tr_uvw, G, dir, bc_tri)) continue;
                    if (Tetrahedron(tr_uvw[0], tr_uvw[1], tr_uvw[2], G).volume() > -1e-10) continue;
                    vec3 lI(0, 0, 0);
                    FOR(lv, 3) lI += bc_tri[lv] * tr_uvw[lv];


                    /*DEBUG*/I = vec3(0, 0, 0); FOR(lv, 3) I += bc_tri[lv] * f.vertex(lv).pos();
                    out_facet = f;
                    if (!f.opposite().active())
                        cond[Volume::Cell(hex, cur_hex).facet(lf)] = { I,ToolBox(tet).facet_geom(f).normal() };
                    break;
                };

                if (out_facet == -1) {
                    Trace::alert("I was not able to cross a tet... maybe due to numerical precision ?");
                    break;
                }
                if (!out_facet.opposite().active()) {
                    break;
                }
                cur_tet = out_facet.opposite().cell();
                if (tf.ap.mid != 0 && GPTransitionFunction(out_facet, uvw).ap.mid != 0) plop("Double transition");
                tf = tf.apply(GPTransitionFunction(out_facet, uvw));
            }

        }


        if (objective_hex != -1) {// link both cur_hex with objective_hex
            //cur_hex
            FOR(di_cur, 2)FOR(dj_cur, 2)FOR(dk_cur, 2) {
                int lv_cur = di_cur + 2 * dj_cur + 4 * dk_cur;
                vec3 pos_v_cur = -vec3(1, 1, 1) + 2 * vec3(di_cur, dj_cur, dk_cur);
                FOR(di_obj, 2)FOR(dj_obj, 2)FOR(dk_obj, 2) {
                    int lv_obj = di_obj + 2 * dj_obj + 4 * dk_obj;
                    vec3 pos_v_obj = -vec3(1, 1, 1) + 2 * vec3(di_obj, dj_obj, dk_obj);
                    pos_v_obj = tf.ap.get_mat() * pos_v_obj;
                    pos_v_obj += 2 * dir;
                    if ((pos_v_obj - pos_v_cur).norm2() == 0)
                        ds.merge(hex.vert(cur_hex, lv_cur), hex.vert(objective_hex, lv_obj));
                }
            }

            plrot[ToolBox(pl).add_segment(pos[cur_hex], pos[objective_hex])] = tf.ap.mid;

        }
        else {// it left the volume on I
            plrot[ToolBox(pl).add_segment(I, pos[cur_hex])] = tf.ap.mid;
        }

    }
    if (verbose > 1) Drop(pl, plrot).apply_wireframe("links");


    ToolBox(hex).merge_vertices(ds);

    hex.connect();
    for (auto f : hex.iter_facets()) if (f.on_boundary())FOR(lv, 4) modified[f.vertex(lv)] = true;

    if (verbose > 1) DropVolume(hex).apply("Connected");

    if (verbose > 1)  show_constraints("constraints");
    return*this;
}







// get projected dimension, position and anisotropy
std::tuple<int, vec3, mat3x3> DualContour::get_fitting_constraint(std::vector<vec3>& pt, std::vector<vec3>& n) {
    mat3x3 AtA;
    vec3 AtB;

    vec3 G;
    FOR(p, pt.size()) G += pt[p];
    G /= double(pt.size());

    auto add_plane = [&](vec3 P, vec3 N) {
        FOR(i, 3)       FOR(j, 3)AtA[i][j] += N[i] * N[j];
        FOR(i, 3)       AtB[i] += P * N * N[i];
        };

    FOR(p, pt.size())   add_plane(pt[p], n[p]);
    FOR(i, 3)           AtA[i][i] += .00001;
    FOR(i, 3)           AtB[i] += .00001 * G[i];

    auto [eigen_val, eigen_vect] = eigendecompose_symmetric(AtA);

    std::tie(eigen_val, eigen_vect) = eigendecompose_symmetric(AtA);
    int dim = 0;
    if (eigen_val[1] < .2 * eigen_val[0]) {
        dim = 2;
        vec3 normal = eigen_vect.col(0);// or column ? idk
        vec3 x = cross(vec3(1, 0, 0), normal);
        if (x.norm2() < .1) x = cross(vec3(0, 1, 0), normal);
        x.normalize();
        vec3 y = cross(normal, x);
        add_plane(G, 10 * x);
        add_plane(G, 10 * y);
        std::tie(eigen_val, eigen_vect) = eigendecompose_symmetric(AtA);
    }
    else
        if (eigen_val[2] < .2 * eigen_val[0]) {
            dim = 1;
            vec3 x = cross(eigen_vect.col(0), eigen_vect.col(1));
            add_plane(G, 10 * x);
            std::tie(eigen_val, eigen_vect) = eigendecompose_symmetric(AtA);
        }
    mat3x3 inv_val;
    FOR(i, 3) inv_val[i][i] = 1. / eigen_val[i];
    mat3x3 eigen_inv = eigen_vect * inv_val * eigen_vect.transpose();

    return { dim,eigen_inv * AtB,eigen_vect };

}



DualContour& DualContour::postpad_propagate_constraints() {
    for (auto f : hex.iter_facets()) if (f.on_boundary() && cond[f].n.norm2() == 0) {
        cond[f] = { Quad3(f).bary_verts(),Quad3(f).normal() };
        FOR(lv, 4) modified[f.vertex(lv)] = true;
    }
    return*this;
}


DualContour& DualContour::show_constraints(std::string name) {
    PolyLine pl;
    double ave = ToolBox(hex).ave_edge_size();
    EdgeAttribute<double> norm(pl);
    for (auto f : hex.iter_facets()) if (f.on_boundary() && cond[f].n.norm2() != 0) {
        vec3 n = .5 * ave * cond[f].n.normalized();
        int e = ToolBox(pl).add_segment(cond[f].pt, cond[f].pt + n);
        norm[e] = cond[f].n.norm();
    }
    if (verbose > 1) Drop(pl, norm).apply_wireframe(name);
    return*this;
}



DualContour& DualContour::pad_for_quad_corner() {
    for (auto v : hex.iter_vertices()) modified[v] = false;

    HexPad pad(hex);
    CellFacetAttribute<bool> to_pad(hex, false);

    for (auto f : hex.iter_facets()) if (f.on_boundary()) {
        vec3 normal = Quad3(f).normal();
        for (auto h : f.iter_halfedges()) {
            vec3 e0 = h.to().pos() - h.from().pos();
            vec3 e1 = h.next().to().pos() - h.to().pos();
            e0.normalize();
            e1.normalize();
            if (cross(e0, e1) * normal<0 || e0 * e1>std::cos(M_PI / 8.)) to_pad[f] = true;
        }
    }

    for (auto c : hex.iter_cells()) {
        bool need_split = false;
        for (auto f : c.iter_facets()) need_split = need_split || (f.on_boundary() && to_pad[f]);
        if (!need_split) continue;
        for (auto f : c.iter_facets()) to_pad[f] = !to_pad[f];
    }
    if (verbose > 1) Drop(hex, to_pad)._skip_value(false).apply("to_pad_deep");

    for (auto f : hex.iter_facets()) if (to_pad[f]) FOR(lv, 4) modified[f.vertex(lv)] = true;

    pad.apply(to_pad, false, false, true);
    postpad_propagate_constraints();
    return*this;
}





DualContour& DualContour::pad_for_hex_edge() {
    //    for(auto v:hex.iter_vertices()) modified[v] = false;

    hex.connect();
    CellFacetAttribute<bool> to_pad(hex, true);
    HexPad pad(hex);

    {
        EdgeGraph eg(hex);
        HexBoundary bound(hex);
        CellAttribute<bool> dont_pad(hex, false);

        bool stable = false;
        while (!stable) {
            stable = true;
            PointAttribute<bool> remove_vert(hex, false);
            // flag edges around faces opposite by only 3 vertices
            for (auto h : bound.quad.iter_halfedges()) {
                if (h.prev().from() == h.opposite().next().to()) {
                    remove_vert[h.from()] = true;
                    remove_vert[h.to()] = true;
                }
            }
            //  lock cells as needed 
            for (auto h : hex.iter_halfedges()) if (remove_vert[h.from()] && !dont_pad[h.cell()]) {
                dont_pad[h.cell()] = true;
                stable = false;
            }
            if (verbose > 1)  Drop(hex, dont_pad)._skip_value(false).apply("dont_pad");


            // find padding 
            for (auto f : hex.iter_facets()) {
                to_pad[f] = true;
                if (dont_pad[f.cell()]) to_pad[f] = false;
                if (!f.on_boundary() && !dont_pad[f.opposite().cell()]) to_pad[f] = false;
            }

            //extract padding
            Quads quad;
            quad.points = hex.points;
            for (auto f : hex.iter_facets()) if (to_pad[f]) {
                auto newf = quad.create_facets(1);
                FOR(lv, 4) quad.vert(newf, lv) = f.vertex(lv);
            }
            quad.connect();
            for (auto v : quad.iter_vertices()) {
                if (!v.halfedge().active()) continue;
                int valence = 0;
                for (auto h : v.iter_halfedges()) valence++;
                auto h = v.halfedge();
                int val_ombrella = 0;
                do { h = h.prev().opposite(); if (!h.active()) break; val_ombrella++; } while (h != v.halfedge());
                if (val_ombrella != valence) {
                    Trace::alert("val_ombrella != valence");
                    remove_vert[v] = true;
                }
            }
            DropSurface(quad).apply("quads to pad");

            // check padding manifoldness
            EdgeAttribute<int> pad_valence(eg, 0);
            for (auto h : hex.iter_halfedges())
                if (to_pad[h.facet()]) pad_valence[eg.edge_from_halfedge(h)]++;
            if (verbose > 1) Drop<PolyLine, EdgeAttribute<int> >(eg, pad_valence).apply_wireframe("pad_valence");

            //  lock cells as needed 
            for (auto e : eg.iter_edges()) if (pad_valence[e] > 1) {
                remove_vert[e.from()] = true;
                remove_vert[e.to()] = true;
            }
            for (auto h : hex.iter_halfedges()) if (remove_vert[h.from()] && !dont_pad[h.cell()]) {
                dont_pad[h.cell()] = true;
                stable = false;
            }
            if (verbose > 1)  Drop(hex, dont_pad)._skip_value(false).apply("dont_pad");
            if (verbose > 1)  Drop(hex, to_pad)._skip_value(false).apply("2pad");
        }
    }


    //int nverts_backup = hex.nverts();
    for (auto f : hex.iter_facets()) if (to_pad[f]) FOR(lv, 4) modified[f.vertex(lv)] = true;
    pad.apply(to_pad, false, false, true);
    postpad_propagate_constraints();
    return*this;
}



void merge_vertices(Hexahedra& hex, DisjointSet& ds, PointAttribute<bool>& modified) {
    Trace::step("merge vertices");
    // merge vertices
    PointAttribute<int> set_id(hex.points);
    int new_nverts = ds.get_sets_id(set_id.ptr->data);
    std::vector<vec4> new_pts(new_nverts);
    FOR(v, hex.nverts()) {
        modified[set_id[v]] = modified[v];
        FOR(d, 3) new_pts[set_id[v]][d] += hex.points[v][d];
        new_pts[set_id[v]][3] += 1.;
    }
    FOR(c, hex.ncells()) FOR(lv, 8) hex.vert(c, lv) = set_id[hex.vert(c, lv)];
    FOR(v, new_nverts) hex.points[v] = vec3(new_pts[v][0], new_pts[v][1], new_pts[v][2]) / new_pts[v][3];
    hex.delete_isolated_vertices();
}


DualContour& DualContour::fix_some_degenerated_cases(bool safe_mode) {
    {
        Trace::step("glue facets sharing 3 vertices");
        DisjointSet ds(hex.nverts());
        {
            HexBoundary bound(hex);
            Quads& quad = bound.quad;
            bool done = false;
            while (!done) {
                done = true;
                for (auto h : quad.iter_halfedges()) {
                    if (ds.root(h.prev().from()) != ds.root(h.opposite().next().to())) continue;
                    int A = h.next().to();
                    int B = h.opposite().prev().from();
                    if (ds.root(A) != ds.root(B)) done = false;
                    ds.merge(A, B);
                    modified[A] = true;
                    modified[B] = true;
                }
            }
        }
        merge_vertices(hex, ds, modified);
        hex.connect();
    }
    if (verbose > 1)  Drop(hex.points, modified)._skip_value(false).apply("in");
    if (safe_mode) return *this;
    {
        Trace::step("complete one ring on valence 4 edges");
        DisjointSet ds(hex.nverts());
        {
            EdgeGraph eg(hex);
            EdgeAttribute<int> val(eg, 0);
            for (auto h : hex.iter_halfedges()) val[eg.edge_from_halfedge(h)]++;

            HexBoundary bound(hex);
            Quads& quad = bound.quad;
            bool done = false;
            while (!done) {
                done = true;
                for (auto h : quad.iter_halfedges()) {
                    if (val[eg.edge_from_halfedge(bound.hex_halfedge(h))] != 4) continue;
                    int A = h.next().to();
                    int B = h.opposite().prev().from();
                    if (ds.root(A) != ds.root(B)) done = false;
                    ds.merge(A, B);
                    modified[A] = true;
                    modified[B] = true;
                    A = h.prev().from();
                    B = h.opposite().next().to();
                    if (ds.root(A) != ds.root(B)) done = false;
                    ds.merge(A, B);
                    modified[A] = true;
                    modified[B] = true;
                }
            }
        }
        merge_vertices(hex, ds, modified);
        hex.connect();
        for (auto f : hex.iter_facets()) if (!f.on_boundary()) cond[f].n = vec3(0, 0, 0);
    }
    if (verbose > 1)  Drop(hex.points, modified)._skip_value(false).apply("in");

    {
        Trace::step("remove hexes that are borderline");
        std::vector<bool> to_kill(hex.ncells(), false);
        auto boundary = [&](Volume::Facet f) {
            if (to_kill[f.cell()]) return false;
            if (!f.opposite().active()) return true;
            return bool(to_kill[f.opposite().cell()]);
            };
        double ave = ToolBox(hex).ave_edge_size();
        while (true) {
            bool need_more_iterations = false;
            for (auto h : hex.iter_halfedges()) {
                auto opp = h.opposite_f();
                auto prev_opp = h.prev().opposite_f();

                // filter h tq h.from() is a convex corner
                if (!boundary(h.facet())) continue;
                if (!boundary(opp.facet())) continue;
                if (!boundary(prev_opp.facet())) continue;


                auto h_data = cond[h.facet()];
                auto opp_data = cond[opp.facet()];
                auto prev_opp_data = cond[prev_opp.facet()];

                // check it the hex is really ambiguous
                const double alpha = M_PI / 4.;
                double normal_threshold = std::pow(std::cos(alpha), 2) + std::pow(1 - std::sin(alpha), 2);
                if ((h_data.pt - opp_data.pt).norm2() > std::pow(.3 * ave, 2)) continue;
                if ((h_data.n - opp_data.n).norm2() > normal_threshold)       continue;
                if ((h_data.n - prev_opp_data.n).norm2() > normal_threshold)  continue;

                Volume::Facet oppf_opp = h.next().next().opposite_f().facet().opposite();
                Volume::Facet oppf_prev_opp = h.next().opposite_f().facet().opposite();
                Volume::Facet oppf_h = opp.next().next().opposite_f().facet().opposite();
                if (oppf_h.active()) cond[oppf_h] = h_data;
                if (oppf_opp.active()) cond[oppf_opp] = opp_data;
                if (oppf_prev_opp.active()) cond[oppf_prev_opp] = prev_opp_data;

                cond[h.facet()].n = vec3(0, 0, 0);
                cond[opp.facet()].n = vec3(0, 0, 0);
                cond[prev_opp.facet()].n = vec3(0, 0, 0);

                to_kill[h.cell()] = true;
                need_more_iterations = true;
            }

            if (!need_more_iterations) break;
        }
        for (auto c : hex.iter_cells()) if (to_kill[c]) FOR(lv, 4) modified[c.vertex(lv)] = true;

        hex.delete_cells(to_kill); hex.delete_isolated_vertices();
        hex.connect();
        if (verbose > 1) { HexBoundary bound(hex); DropSurface(bound.quad).apply("boundary"); }
    }


    show_constraints("constraints After fix");
    return*this;
}


using namespace UM::Linear;

DualContour& DualContour::compute_charts(HexBoundary& bound,FacetAttribute<int>& chart,int& nb_charts){
    Quads& quads = bound.quad;
    //Hexahedra& hex = bound.hex;

    CornerAttribute<int> nhex(quads,0);
    {
        EdgeGraph eg(hex);
        EdgeAttribute<int> edge_valence(eg, 0);
        for (auto h : hex.iter_halfedges())     edge_valence[eg.edge_from_halfedge(h)]++;
        for (auto h : quads.iter_halfedges())   nhex[h] = edge_valence[eg.edge_from_halfedge(bound.hex_halfedge(h))];//1;
    }

    PointAttribute<int> val(quads,0);
    for (auto h : quads.iter_halfedges())  val[h.from()]++;

    CornerAttribute<bool> cut(quads,false);
    for (auto h : quads.iter_halfedges())  cut[h] = 2!=nhex[h] ;

    for (auto seed : quads.iter_halfedges()){
        if (!cut[seed]) continue;
        if (!cut[seed.prev()]) continue;
        if (!cut[seed.opposite().next()]) continue;
        std::vector<int> cut_stack;
        Surface::Halfedge h = seed;
        bool valid_cut = true;
        while (valid_cut){
            cut_stack.push_back(h);
            if (!cut[h])                                                {valid_cut = false; break;}
            if (cut[h.next()] != cut[h.opposite().prev()])              {valid_cut = false; break;}

            if (cut[h.next()] && cut[h.opposite().prev()])              break; // success
            
            h = h.next().opposite().next();
            if (val[h.from()]!=4)                                       {valid_cut = false; break;}

        } 
        if (!valid_cut) continue;
        valid_cut = false;
        for(auto hid: cut_stack) {
            Surface::Halfedge lh(quads,hid);
            if (cond[bound.hex_facet(lh.facet())].n*cond[bound.hex_facet(lh.opposite().facet())].n>std::cos(M_PI/4.))
                valid_cut = true;
        }
        if (!valid_cut) continue;
        for(auto hid: cut_stack) {cut[hid]=false;cut[Surface::Halfedge(quads,hid).opposite()]=false;}
    }
    
    // render
    if (verbose > 1)  Drop(quads,nhex)._wireframe(true).apply_half_edge("nhex");
    if (verbose > 1) Drop(quads,cut)._skip_value(false)._wireframe(true).apply_half_edge("cut");

    auto corner_angle = [&](Surface::Halfedge h){
        vec3 v0 = Segment3(h).vector().normalized();
         vec3 v1 = -Segment3(h.prev()).vector().normalized();
         return std::atan2(cross(v0,v1).norm(),v0*v1);
    };

    // cut concave
    while (true){
        std::vector<int> best_cut;
        for(auto seed:quads.iter_halfedges()){
            if (cut[seed]) continue;
            if (!cut[seed.prev()]) continue;
            if (cut[seed.opposite().next()]) continue;            

            if (corner_angle(seed) + corner_angle(seed.opposite().next()) + corner_angle(seed.opposite().next().opposite().next())<1.2*M_PI) continue;

            for(int hid:{seed,seed.opposite().next()}){
                std::vector<int> cur_cut;
                Surface::Halfedge iter(quads, hid);
                bool valid_cut = true;
                while (true) {
                    cur_cut.push_back(iter);
                    if (cut[iter]) { valid_cut = false; break; }
                    if (cut[iter.next()] || cut[iter.opposite().prev()]) break; 
                    iter = iter.next().opposite().next();
                    if (val[iter.from()] != 4) { valid_cut = false; break; }
                }
                //plop(best_cut.size());
                if (!valid_cut) continue;
                if (best_cut.empty() || cur_cut.size()<best_cut.size())
                    best_cut = cur_cut;            
            }
        }
        //plop(best_cut.size());
        if (best_cut.empty()) break;
        for(auto hid: best_cut) {cut[hid]=true;cut[Surface::Halfedge(quads,hid).opposite()]=true;}
    }


    if (verbose > 1) Drop(quads,cut)._skip_value(false)._wireframe(true).apply_half_edge("cuttmeshified");


        // very naive implementation: disjointset with frontiers detected by dihedral angle
    { 
        DisjointSet ds(quads.nfacets());
        for (auto h : quads.iter_halfedges())
            if (cond[bound.hex_facet(h.facet())].n * cond[bound.hex_facet(h.opposite().facet())].n  > std::cos(M_PI / 6.))
                ds.merge(h.facet(), h.opposite().facet());
        nb_charts = ds.get_sets_id(chart.ptr->data);
    }

    return*this;
}



bool DualContour::try_to_pad() {
    for (auto v : hex.iter_vertices()) modified[v] = false;

    HexPad pad(hex);
    CellFacetAttribute<bool> to_pad(hex, false);

    {

        HexBoundary bound(hex);
        Quads& quads = bound.quad;
        FacetAttribute<int> chart(quads);
        int nb_charts;


        // init charts
        compute_charts(bound,chart,nb_charts);
        if (verbose > 1) Drop(quads, chart).apply("charts");

        //return false;





        // init edge graph attributes
        EdgeGraph eg(hex);
        EdgeAttribute<int> edge_valence(eg, 0);
        EdgeAttribute<bool> edge_on_border(eg, false);
        EdgeAttribute<bool> inside_chart(eg, false);
        for (auto h : hex.iter_halfedges()) {
            auto e = eg.edge_from_halfedge(h);
            edge_valence[e]++;
            if (h.facet().on_boundary()) edge_on_border[e] = true;
        }
        for (auto h : quads.iter_halfedges()) {
            if (chart[h.facet()] == chart[h.opposite().facet()])
                inside_chart[eg.edge_from_halfedge(bound.hex_halfedge(h))] = true;
        }

        std::vector<bool> need_pad(nb_charts, false);
        for (auto h : quads.iter_halfedges()) {
            if (chart[h.facet()] != chart[h.opposite().facet()]) continue;
            if (edge_valence[eg.edge_from_halfedge(bound.hex_halfedge(h))] == 1)
                need_pad[chart[h.facet()]] = true;
        }


        FOR(cur_chart, nb_charts) {
            plop(cur_chart);
            if (!need_pad[cur_chart]) continue;
            for (auto f : hex.iter_facets()) to_pad[f] = false;
            for (auto f : quads.iter_facets()) if (chart[f] == cur_chart) to_pad[bound.hex_facet(f)] = true;

            // expand charts inside
#if 1 //simply cut straight
            bool done = false;
            while (!done) {
                done = true;
                for (auto f : hex.iter_facets()) if (to_pad[f]) {
                    for (auto f_cir : f.iter_halfedges()) {
                        auto e = eg.edge_from_halfedge(f_cir);
                        if (edge_valence[e] != 4 && !edge_on_border[e]) continue;
                        if (inside_chart[e]) continue;
                        auto next = f_cir.opposite_f().opposite_c();
                        if (!next.active()) continue;
                        if (!to_pad[next.opposite_f().facet()]) {
                            to_pad[next.opposite_f().facet()] = true;
                            done = false;
                        }
                    }
                }
            }
#endif


            // check flagging validity
            bool flag_is_valid = true;
            EdgeAttribute<int> npads(eg, 0);
            for (auto h : hex.iter_halfedges())
                if (to_pad[h.facet()])
                    npads[eg.edge_from_halfedge(h)]++;
            for (auto f : quads.iter_facets())
                if (to_pad[bound.hex_facet(f)])
                    for (auto h : f.iter_halfedges()) {
                        auto e = eg.edge_from_halfedge(bound.hex_halfedge(h));
                        if (npads[e] + npads[e.opposite()] == 1 && edge_valence[e] != 1)
                            flag_is_valid = false;
                    }
            if (!flag_is_valid) continue;

            if (verbose > 1)  Drop(hex, to_pad)._skip_value(false).apply("pad_test");
            goto gogogo;
        }
    }
    return false;

gogogo:
    for (auto f : hex.iter_facets()) if (to_pad[f]) FOR(lv, 4) modified[f.vertex(lv)] = true;
    pad.apply(to_pad, false, false, true);
    postpad_propagate_constraints();
    return true;
}
























#include <fullhex/untangle_tetrahedra.h>

using namespace UM::Linear;


UntanglerHexTan::UntanglerHexTan(Hexahedra& m) :m(m) {}

    void UntanglerHexTan::apply(DualContour& dc, double fit_scale) {
        // create the active region
        PointAttribute<int>  dist_to_modified(m, 10000);
        int n_ring_size = 2;
        {
            for (auto v : m.iter_vertices()) if (dc.modified[v]) dist_to_modified[v] = 0;

            FOR(i, n_ring_size) for (auto h : m.iter_halfedges())
                if (dist_to_modified[h.to()] > dist_to_modified[h.from()] + 1)
                    dist_to_modified[h.to()] = dist_to_modified[h.from()] + 1;
        }


        Tetrahedra tets;
        tets.points = m.points;
        for (auto c : m.iter_cells()) FOR(lv, 8) {
            FOR(lh, c.nhalfedges()) {
                auto h = c.halfedge(lh);
                if (h.from() != c.vertex(lv)) continue;

                int verts[4] = { h.from(),h.to(),h.opposite_f().next().to(),h.prev().from() };
                bool completely_locked = true;
                FOR(lv0, 4)  completely_locked = completely_locked && dist_to_modified[verts[lv0]] > n_ring_size;

                if (completely_locked) break;

                int nt = tets.create_cells(1);
                FOR(lv1, 4)  tets.vert(nt, lv1) = verts[lv1];

                break;
            }
        }
        //plop(tets.ncells());

        //DropVolume(tets).apply("untangletets");
        tets.connect();
        UntanglerTet untangle(tets);
        for (auto v : m.iter_vertices()) untangle.lock[v] = dist_to_modified[v] > n_ring_size;

        double a = ToolBox(untangle.m).ave_edge_size();
        for (auto t : tets.iter_cells()) {
            untangle.volume[t] = a * a * a / 6.;//vol_ave;
            vec3 tet[4];
            tet[0] = a * vec3{ 0,0,0 };
            tet[1] = a * vec3{ 1,0,0 };
            tet[2] = a * vec3{ 0,1,0 };
            tet[3] = a * vec3{ 0,0,1 };
            untangle.reference[t] = mat<4, 3>{ { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} } }*mat<3, 3>{{tet[1] - tet[0], tet[2] - tet[0], tet[3] - tet[0]}}.invert_transpose();
        }

        PointAttribute< std::vector<vec3> > v_pt(m);
        PointAttribute< std::vector<vec3> > v_n(m);

        for (auto f : dc.hex.iter_facets()) if (f.on_boundary())
            for (auto h : f.iter_halfedges()) {
                v_pt[h.from()].push_back(dc.cond[f].pt);
                v_n[h.from()].push_back(dc.cond[f].n);
            }

        PointAttribute<int> pdim(m, -1);
        PointAttribute<vec3> firsteigen(m, vec3(0, 0, 0));
        for (auto v : m.iter_vertices()) {
            if (v_n[v].empty()) continue;
            auto [dim, P, eigen_vect] = dc.get_fitting_constraint(v_pt[v], v_n[v]);
            pdim[v] = dim;
            firsteigen[v] = eigen_vect.col(2);
            FOR(ax, 3) {
                LinExpr line = -P * eigen_vect.col(ax);
                FOR(d, 3) line += eigen_vect.col(ax)[d] * X(3 * v + d);
                if (ax >= dim)     untangle.ls_matrix.push_back(fit_scale * 10 * line);
                else                untangle.ls_matrix.push_back(fit_scale * line);
            }
        }
        untangle.apply(false);
    }
