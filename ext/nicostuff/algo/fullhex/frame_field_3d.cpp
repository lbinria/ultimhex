
#include <drop_attribute.h>
#include <framework/trace.h>
#include <fullhex/frame_field_3d.h>
#include <fullhex/gp_basic.h>


FF3D::FF3D(Tetrahedra& m, CellAttribute<mat3x3>& J) : m(m), ap(m), J(J) {}

mat<3, 3>& FF3D::operator[](int i) { return J[i]; }


AxisPermutation FF3D::permutation_around_edge(Volume::Halfedge h) {
    auto cir = h;
    AxisPermutation perm(0);
    do {
        //perm = ap[cir.opposite_f().facet()] * perm;
        perm = perm * ap[cir.opposite_f().facet()];
        cir = cir.opposite_f().opposite_c();
        if (cir == -1) return false;
    } while (cir != h);
    return perm;
}

bool FF3D::edge_is_singular(Volume::Halfedge h) {
    return !permutation_around_edge(h).is_identity();
}

int FF3D::singular_edge_stable_coordinate(Volume::Halfedge h) {
    return permutation_around_edge(h).stable_axis();
}

bool FF3D::trace_streamline(int root_cell, int branch, vec3 P, std::vector<int>& cells, std::vector<vec3>& pts, int max_length) {
    cells.push_back(root_cell);
    pts.push_back(P);
    vec3 dir_uvw = label_to_normal[branch];
    FOR(it, max_length) {
        Volume::Cell c(m, cells.back());
        vec3 dir = J[c].transpose() * dir_uvw;
        P = pts.back();
        for (auto f : c.iter_facets()) {
            vec3 bc;
            Triangle3 tr = Triangle3(f);
            if (Intersect::triangle_line(tr, P, dir, bc)) {
                auto opp = f.opposite();
                if (opp.active() && cells.size() > 2 && opp.cell() == cells[cells.size() - 2]) continue;
                pts.push_back(bc[0] * tr[0] + bc[1] * tr[1] + bc[2] * tr[2]);
                if (!opp.active()) break;
                cells.push_back(opp.cell());
                AxisPermutation l_ap;
                l_ap = permute_Jj_to_approx_Ji(J[opp.cell()], J[c]);
                dir_uvw = l_ap.get_mat() * dir_uvw;
                break;
            }
        }
        if (c == cells.back() || !Volume::Cell(m, cells.back()).active())  break;
    }
    return cells.size() < max_length;
}


void FF3D::show_streamline(std::string name, int nb_shots) {
    PolyLine pl;
    EdgeAttribute<int> path(pl, 0);
    int num_path = 0;
    FOR(shot, nb_shots) {
        Volume::Cell root_cell(m, rand() % m.ncells());

        int branch_seed = rand() % 3;
        FOR(d, 2) {
            //int max_length = 100;
            int branch = d + 2 * branch_seed;
            vec3 P = Tetrahedron(root_cell).bary_verts();
            std::vector<int> cells;
            std::vector<vec3> pts;
            trace_streamline(root_cell, branch, P, cells, pts, 1000);
            FOR(i, pts.size() - 1) path[ToolBox(pl).add_segment(pts[i], pts[i + 1])] = num_path;
        }
        num_path++;
    }
    Drop(pl, path).apply_wireframe(name);
}
void FF3D::show_singularity_graph(std::string name, bool with_neig_cubes) {
    EdgeGraph eg(m);
    EdgeAttribute<bool> singu(eg);
    for (auto e : eg.iter_edges()) singu[e] = edge_is_singular(eg.halfedge_from_edge(e));
    Drop<PolyLine, EdgeAttribute<bool>  >(eg, singu)._skip_value(false).apply_wireframe(name);


    if (with_neig_cubes) {
        double ave = ToolBox(m).ave_edge_size();
        Hexahedra outm;
        CellFacetAttribute<int> flag(outm, 0);

        int nc = 0;
        for (auto e : eg.iter_edges()) if (singu[e])
            for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge())
                nc++;

        outm.create_cells(nc);
        outm.points.create_points(8 * nc);
        FOR(c, nc)	FOR(f, 6) flag[outm.facet(c, f)] = f;
        FOR(c, nc)	FOR(v, 8) outm.vert(c, v) = 8 * c + v;


        int outc = 0;
        for (auto e : eg.iter_edges()) if (singu[e])
            for (auto cir : eg.halfedge_from_edge(e).iter_CCW_around_edge()) {
                auto c = cir.cell();
                FOR(di, 2)FOR(dj, 2)FOR(dk, 2)
                    outm.points[8 * outc + 4 * di + 2 * dj + dk] = Tetrahedron(c).bary_verts() + 0.25 * ave * J[c].transpose()/*J is a rot => transpose == fast_invert*/ * vec3(2 * di - 1, 2 * dj - 1, 2 * dk - 1);

                outc++;
            }

        if (outm.nfacets() > 0)
            Drop(outm, flag).apply(name + "_cubes");
    }


}
void FF3D::show_cubes(std::string name, double scale, bool with_flag) {
    double ave = ToolBox(m).ave_edge_size();
    Hexahedra outm;
    outm.create_cells(m.ncells());
    outm.points.create_points(8 * m.ncells());
    CellFacetAttribute<int> flag(outm, 0);
    FOR(f, outm.nfacets())flag[f] = f % 6;
    for (auto c : m.iter_cells()) {
        FOR(v, 8) outm.vert(c, v) = 8 * c + v;
        FOR(di, 2)FOR(dj, 2)FOR(dk, 2) outm.points[8 * c + di + 2 * dj + 4 * dk] = Tetrahedron(c).bary_verts() + scale * .15 * ave
            * J[c].transpose() /*J is a rot => transpose == fast_invert*/ * vec3(2 * di - 1, 2 * dj - 1, 2 * dk - 1);
    }
    if (with_flag)	Drop(outm, flag).apply(name + "_cubes");
    else			DropVolume(outm).apply(name + "_cubes");
}
void FF3D::show(std::string prefix) {
    if (!Trace::drop_mesh_is_active) return;
    //show_cubes(prefix+ "cube");
    show_singularity_graph(prefix + "singu", true);
    show_streamline(prefix + "stream");
}



void FF3D::init_with_constant_frame() {
    CellFacetAttribute<double> quality(m, 1);
    std::vector<vec3> bunch_of_vectors;
    for (auto f : m.iter_facets()) if (!f.opposite().active()) bunch_of_vectors.push_back(Triangle3(f).normal());
    mat<3, 3> M = Frame::representative_frame(bunch_of_vectors);
    for (int c : m.iter_cells()) J[c] = M;
}
void FF3D::init_from_uvw(CellCornerAttribute<vec3>& U, bool smooth) {
    map_to_jacobian(m, U, J);
    CellAttribute<mat3x3> accum(m);
    std::copy(J.ptr->data.begin(), J.ptr->data.end(), accum.ptr->data.begin());
    for (auto f : m.iter_facets()) if (!f.on_boundary()) accum[f.cell()] += J[f.opposite().cell()];
    std::copy(accum.ptr->data.begin(), accum.ptr->data.end(), J.ptr->data.begin());
    for (auto c : m.iter_cells()) J[c] = closest_rotation(J[c]);
}
void FF3D::update_axis_permutation() {
    for (auto f : m.iter_facets()) {
        auto opp = f.opposite();
        if (opp == -1)	ap[f].mid = 0;
        if (opp < f) continue;
        ap[f] = permute_Jj_to_approx_Ji(J[f.cell()], J[opp.cell()]);
        ap[opp] = ap[f].inverse();
    }
}



using namespace UM::Linear;

void FF3D::compute_with_fibers(CellAttribute<int>& fiber, CellFacetAttribute<int>& constraint_type, CellAttribute<bool>& locked, double data_fitting_w, int nb_pass) {
    //um_assert(nb_pass ==1); // multipass is likely to be better, BUT: 1- slower, 2- not tested yet

    CellAttribute<int> rep_fiber(m, -1);
    for (auto c : m.iter_cells()) if (fiber[c] != -1) rep_fiber[fiber[c]] = c;

    CellAttribute<double> angle(m, 0);



    FOR(iter, nb_pass) {
        ConstrainedLeastSquares cls(2 * m.ncells());


        // constant along fiber
        for (auto c : m.iter_cells()) if (!locked[c]) FOR(d, 2) cls.add_to_constraints(X(2 * rep_fiber[fiber[c]] + d) - X(2 * c + d));

        // no rotation on locked cells
        for (auto c : m.iter_cells()) if (locked[c]) FOR(d, 2) cls.add_to_constraints(X(2 * c + d) - (1 - d));

        // smoothing term
        for (auto f : m.iter_facets()) {
            if (!f.opposite().active()) continue;
            FOR(l, 2) cls.add_to_energy(X(2 * f.opposite().cell() + l) - X(2 * f.cell() + l));
        }

        //data fitting
        for (auto c : m.iter_cells()) if (iter > 0 && fiber[c] != -1) {
            double fitw = .1;
            //angle[c] += rand_range(-.1,.1);
            cls.add_to_energy(fitw * (X(2 * c + 0) - std::cos(4. * angle[c])));
            cls.add_to_energy(fitw * (X(2 * c + 1) - std::sin(4. * angle[c])));
        }

        // boundary fitting
        if (data_fitting_w > 0) for (auto f : m.iter_facets()) {
            if (f.opposite().active()) continue;
            vec3 n = Triangle3(f).normal();
            auto c = f.cell();

            vec3 n_B = J[c] * n;

            if (constraint_type[f] != 5) continue; //(boundary is free or perp to z)

            double a = 4. * std::atan2(n_B[1], n_B[0]);
            cls.add_to_energy(data_fitting_w * (X(2 * c) - std::cos(a)));
            cls.add_to_energy(data_fitting_w * (X(2 * c + 1) - std::sin(a)));
        }

        cls.solve();
        for (auto c : m.iter_cells())  angle[c] = std::atan2(cls.value(2 * c + 1), cls.value(2 * c)) / 4.;
    }
    for (auto c : m.iter_cells()) J[c] = rotz(-angle[c]) * J[c];


    update_axis_permutation();
}





bool FF3D::optimize_topology() {
    Trace::Section sect("optimize_topology");
    update_axis_permutation();
    show_singularity_graph("inSG", false);

    DropVolume(m).apply("geom");
    double ave = ToolBox(m).ave_edge_size();
    EdgeGraph eg(m);
    bool done = false;
    while (!done) {
        done = true;
        for (auto v : eg.iter_vertices()) {

            EdgeGraph::Edge e0(eg, -1);
            EdgeGraph::Edge e1(eg, -1);
            {

                for (auto e : v.iter_edges()) {
                    if (!edge_is_singular(eg.halfedge_from_edge(e))) continue;//singu_edge[e]==-1) continue;
                    if (!e0.active()) e0 = e;
                    else if (!e1.active()) e1 = e;
                    else {
                        std::cerr << "FF failed because high (>2) valence node detected in singu graph\n";
                        return false;
                        e1 = -1; //break;
                        show("Fail");
                        throw (std::runtime_error("Fail"));
                        um_assert(!"high (>2) valence node detected in singu graph");
                    }
                }
                if (!e1.active()) continue;
            }

            std::vector<int> path;
            {
                std::map<int, int> prev;
                std::map<int, double> dist;

                auto cmp = [&](int left, int right) { return dist[left] > dist[right]; };
                std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);

                for (auto cir : eg.halfedge_from_edge(e0).iter_CCW_around_edge()) {
                    dist[cir] = 0;
                    prev[cir] = -1;
                    queue.push(cir);
                }

                while (!queue.empty()) {
                    Volume::Halfedge h(m, queue.top());
                    queue.pop();
                    if (eg.edge_from_halfedge(h) == e1) {
                        //plop("youpi ! ");
                        while (true) {
                            h = prev[h];
                            //plop(h);
                            if (h == -1) break;
                            path.push_back(h.next());
                        }
                        break;
                    }

                    auto next = h.prev().opposite_f();
                    double new_path_dist = dist[h] + (h.next().to().pos() - h.to().pos()).norm();
                    if (dist.find(next) == dist.end()) dist[next] = 1e20;
                    if (dist[next] > new_path_dist) {
                        for (auto cir : next.iter_CCW_around_edge()) {

                            dist[cir] = new_path_dist;
                            prev[cir] = h;
                            queue.push(cir);
                        }
                    }
                }
            }

            if (path.empty()) continue;
            double path_length = 0;
            FOR(i, path.size()) path_length += (Volume::Halfedge(m, path[i]).from().pos() - Volume::Halfedge(m, path[i]).to().pos()).norm();
            double actual_length = (e0.from().pos() - e0.to().pos()).norm() + (e1.from().pos() - e1.to().pos()).norm();

            if (actual_length < path_length + .01 * ave) continue;

            done = false;


            FOR(l_h_id, path.size()) {
                Volume::Halfedge h(m, path[l_h_id]);

                FOR(i, 24) {
                    ap[h.facet()].mid = i;
                    ap[h.facet().opposite()] = ap[h.facet()].inverse();
                    if (!edge_is_singular(h.next()))  break;
                }
            }

        }
    }
    show_singularity_graph("outSG", false);
    return true;
}





void FF3D::smooth_geom_only() {
    FOR(it, 10)for (auto c : m.iter_cells()) {
        mat3x3 newJ;
        bool on_border = false;
        FOR(lf, 4) {
            auto f = c.facet(lf);
            auto oppf = f.opposite();
            if (!oppf.active()) {
                on_border = true;
                continue;
            }
            mat<3, 3> rot = ap[f].get_mat() * J[oppf.cell()];
            newJ = newJ + rot;
        }
        if (!on_border)
            J[c] = closest_rotation(newJ);
    }
    update_axis_permutation();
    show();
}



bool FF3D::apply_with_fibers(CellCornerAttribute<vec3>& U, CellAttribute<int>& fibers, CellFacetAttribute<int>& constraint_type, double data_fitting_w) {
    //ff.inflate_ambiguous_boundary_with_fiber(out.U, out.fibers);
//return;
    init_from_uvw(U);
    //update_axis_permutation();
    //show("FF_in");


    CellAttribute<bool> locked(m, false);
    {// lock enough cells to be sure that no singu will touch them
        //for (auto c : m.iter_cells()) locked[c] = (fiber[c] == -1);
        for (auto f : m.iter_facets()) if (constraint_type[f] == -2) locked[f.cell()] = true;

        PointAttribute<bool> vlock(m, false);
        for (auto c : m.iter_cells()) if (locked[c]) FOR(lv, 4) vlock[c.vertex(lv)] = true;
        for (auto c : m.iter_cells()) FOR(lv, 4) if (vlock[c.vertex(lv)])   locked[c] = true;
    }
    Drop(m, locked)._skip_value(false).apply("FFlocked");



    compute_with_fibers(fibers, constraint_type, locked, data_fitting_w, 3);
    //compute_with_fibers(fibers, constraint_type, 1, 3);
    show("FF_geopt");
    TRACE_OFF(optimize_topology());
    TRACE_OFF(smooth_geom_only());
    show("FF_");
    return true;

}


void FF3D::compute_with_fibers_very_smooth(CellAttribute<int>& fiber, CellFacetAttribute<int>& constraint_type) {

    CellAttribute<int> rep_fiber(m, -1);
    for (auto c : m.iter_cells()) if (fiber[c] != -1) rep_fiber[fiber[c]] = c;



    CellAttribute<mat3x3> J_ref(m);
    for (auto c : m.iter_cells()) J_ref[c] = J[c];

    CellAttribute<double> angle(m, 0);



    FOR(iter, 2) {

        CellAttribute<bool> locked(m, false);
        if (iter == 1) {
            // free region close to tangled cells
            PointAttribute<int> pb_dist(m, 10000);
            EdgeGraph eg(m);
            for (auto e : eg.iter_edges())
                if (edge_is_singular(eg.halfedge_from_edge(e))) {
                    pb_dist[e.from()] = 0;
                    pb_dist[e.to()] = 0;
                }
            FOR(it, 10)for (auto e : eg.iter_edges())if (pb_dist[e.from()] > 1 + pb_dist[e.to()])
                pb_dist[e.from()] = 1 + pb_dist[e.to()];
            Drop(m.points, pb_dist).apply("pbdist");
            for (auto f : m.iter_facets()) FOR(lv, 3) if (pb_dist[f.vertex(lv)] > 4) locked[f.cell()] = true;
            Drop(m, locked).apply("lock");
        }

        ConstrainedLeastSquares cls(2 * m.ncells());
        // constant along fiber
        for (auto c : m.iter_cells()) if (fiber[c] >= 0) FOR(d, 2) cls.add_to_constraints(X(2 * rep_fiber[fiber[c]] + d) - X(2 * c + d));


        if (iter == 1) for (auto c : m.iter_cells()) if (locked[c] && (fiber[c] == -1 || rep_fiber[fiber[c]] == c)) { // no rotation on locked cells
            cls.add_to_constraints(X(2 * c + 0) - std::cos(4. * angle[c]));
            cls.add_to_constraints(X(2 * c + 1) - std::sin(4. * angle[c]));
        }

        // smoothing term
        for (auto f : m.iter_facets()) if (!f.on_boundary()) FOR(l, 2)
            cls.add_to_energy(X(2 * f.opposite().cell() + l) - X(2 * f.cell() + l));




        if (iter == 0) for (auto f : m.iter_facets()) {// boundary fitting
            if (f.opposite().active()) continue;
            vec3 n = Triangle3(f).normal();
            auto c = f.cell();

            vec3 n_B = J[c] * n;

            if (constraint_type[f] != 5) continue; //(boundary is free or perp to z)

            //plop(f);
            double a = 4. * std::atan2(n_B[1], n_B[0]);
            cls.add_to_energy(100 * (X(2 * c) - std::cos(a)));
            cls.add_to_energy(100 * (X(2 * c + 1) - std::sin(a)));
        }



        //if (iter == 1) return;
        cls.solve();
        for (auto c : m.iter_cells())  angle[c] = std::atan2(cls.value(2 * c + 1), cls.value(2 * c)) / 4.;
        for (auto c : m.iter_cells()) J[c] = rotz(-angle[c]) * J_ref[c];
        update_axis_permutation();
        show("OPT1");
        show_cubes("cub");
    }



}


bool FF3D::apply_with_fibers_very_smooth(CellCornerAttribute<vec3>& U, CellAttribute<int>& fibers, CellFacetAttribute<int>& constraint_type) {
    init_from_uvw(U);

    compute_with_fibers_very_smooth(fibers, constraint_type);
    //show("FF_");
    //return true;

    TRACE_OFF(optimize_topology());
    TRACE_OFF(smooth_geom_only());

    return true;
}




void FF3D::apply_ls(CellFacetAttribute<bool>& lock_normal, double normal_coeff,CellAttribute<bool>* lock_cells) {
    PointAttribute<SH4> sh(m);
    LeastSquares ls(11 * m.ncells());
    m.connect();

    for (auto f:m.iter_facets()) if (lock_normal[f]) 
        J[f.cell()] = Triangle3(f).tangent_basis().transpose();

    if (lock_cells!=NULL) FOR(c, m.ncells()) if ((*lock_cells)[c]){
        SH4 l_sh;
        l_sh[4] = std::sqrt(7. / 12.);
        l_sh[8] = std::sqrt(5. / 12.);
        vec3 xyz = mat_to_euler(normalize_columns(J[c].transpose()));
        l_sh.euler_rot(xyz);
        FOR(i, 9) ls.fix(11 * c + i, l_sh[i]);
    }

    Trace::step("Construct matrix");
    for (auto f : m.iter_facets()) if (!f.on_boundary())
        FOR(i, 9) ls.add_to_energy(X(11 * f.cell() + i) - X(11 * f.opposite().cell() + i));

    // locked normal on facets enforced by barrier equations
    FOR(c, m.ncells()) {
        if (!lock_normal[m.facet(c, 0)] && !lock_normal[m.facet(c, 1)] && !lock_normal[m.facet(c, 2)] && !lock_normal[m.facet(c, 3)]) continue;
        SH4 sh0, sh4, sh8;
        sh4[4] = std::sqrt(7. / 12.);
        sh0[0] = std::sqrt(5. / 12.);
        sh8[8] = std::sqrt(5. / 12.);
        vec3 xyz = mat_to_euler(normalize_columns(J[c]));
        sh4.euler_rot(xyz);
        sh0.euler_rot(xyz);
        sh8.euler_rot(xyz);
        FOR(i, 9) ls.add_to_energy(normal_coeff * (X(11 * c + i) + sh0[i] * X(11 * c + 9) + sh8[i] * X(11 * c + 10) - sh4[i]));
    }

    Trace::step("Solve SH");
    ls.solve();

    Trace::step("Project SH to B");
    // convert spherical harmonic coefficients to a rotation
    {
//#pragma omp parallel for
        FOR(c, m.ncells()) {
            SH4 fv;
            FOR(i, 9) fv[i] = ls.value(c * 11 + i);
            //FOR(i, 9) std::cerr<<fv[i]<<"  ";std::cerr<<std::endl;
            J[c] = fv.project_mat(1e-3, 1e-5, nullptr).transpose();
        }
    }
     update_axis_permutation();
}
