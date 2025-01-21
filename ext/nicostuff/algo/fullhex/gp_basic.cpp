#include <fullhex/gp_basic.h>
#include <framework/trace.h>

#include <toolbox.h>
#include <drop_attribute.h>



SparseVector express_with_free_variables(ConstrainedLeastSquares& cls, const LinExpr& le) {
    SparseVector v = le;
    v.compact();
    if (!v.empty() && v.front().index < 0)
        v.front().index = cls.rb.C.nrows() - 1; // N.B. it is not sorted anymore
    cls.rb.leading_to_free(v);
    return v;
}

// usage:  auto [already_satisfied,impossible] = constraint_status(cls,le);
std::tuple<bool, bool> constraint_status(ConstrainedLeastSquares& cls, const LinExpr& le) {
    SparseVector v = express_with_free_variables(cls, le);
    if (v.empty())
        return { true,false };
    if (v.size() == 1 && v[0].index == cls.rb.C.nrows() - 1)
        return{ false,true };
    return{ false,false };
}

// jacobian: chaque col est le gradiant d'une coord
mat3x3 uvw_to_jacobian(Volume::Cell c, CellCornerAttribute<vec3>& uvw) {
    mat3x3 J;
    Tetrahedron tet = c;
    mat<3, 4> grd = tet.grad_operator();
    FOR(d, 3) {
        vec4 scal;
        FOR(lv, 4) scal[lv] = uvw[c.corner(lv)][d];
        J[d] = grd * scal;
    }
    return J;
}

mat<4, 3> jacobian_to_uvw(Volume::Cell c, CellAttribute<mat3x3>& J) {
    mat<4, 3> uvw;
    FOR(lv, 4) uvw[lv] = J[c] * (c.vertex(lv).pos() - c.vertex(0).pos());
    return uvw;
}

vec3 uvw(Volume::Corner c, CellAttribute<mat3x3>& J, CellAttribute<vec3>& T) {
    return  J[c.cell()] * (c.vertex().pos() - c.cell().vertex(0).pos()) + T[c.cell()];
}

mat3x3 closest_rotation(mat3x3 J) {
    mat3x3 Jt = J.transpose();
    FOR(d, 3) Jt[d].normalize();
    mat3x3 ortho[3];
    FOR(ld, 3) {
        mat3x3& M = ortho[ld];
        M[0] = Jt[ld];
        M[1] = Jt[(ld + 1) % 3];
        M[2] = Jt[(ld + 2) % 3];
        M[2] = cross(M[0], M[1]).normalized();
        M[1] = cross(M[2], M[0]).normalized();
    }
    FOR(i, 3) Jt[i] = vec3(0, 0, 0);
    FOR(ld, 3) FOR(d, 3) Jt[d] += ortho[ld][(3 - ld + d) % 3];
    FOR(d, 3) Jt[d].normalize();
    return Jt.transpose();
}


AxisPermutation  permute_Jj_to_approx_Ji(mat<3, 3> Ji, mat<3, 3> Jj) {
    int best_i = 0;
    double min_diff = 1e20;
    FOR(i, 24) {
        mat<3, 3> m = AxisPermutation(i).get_mat() * Jj - Ji;
        double diff = m.norm();
        if (diff < min_diff) {
            best_i = i;
            min_diff = diff;
        }
    }
    return AxisPermutation(best_i);
}




void map_to_jacobian(Tetrahedra& m, CellCornerAttribute<vec3>& U_in, CellAttribute<mat3x3>& J_out) {
    for (auto c : m.iter_cells()) {
        Tetrahedron tet = Tetrahedron(c);
        mat<3, 4> grd = tet.grad_operator();
        FOR(d, 3) {
            vec4 scal;
            FOR(lv, 4) scal[lv] = U_in[c.corner(lv)][d];
            J_out[c][d] = grd * scal;
        }
    }

}
void map_to_jacobian(Tetrahedra& m, CellCornerAttribute<vec3>& U_in, CellAttribute<mat3x3>& J_out, CellAttribute<vec3>& T_out) {
    map_to_jacobian(m, U_in, J_out);
    for (auto c : m.iter_cells()) T_out[c] = U_in[c.corner(0)];
}

void jacobian_to_map(Tetrahedra& m, CellAttribute<mat3x3>& J_in, CellAttribute<vec3>& T_in, CellCornerAttribute<vec3>& U_out) {
    for (auto c : m.iter_corners()) U_out[c] = uvw(c, J_in, T_in);
}






GPTransitionFunction::GPTransitionFunction() :t(0, 0, 0) { ap.mid = 0; };


// the transition function express the f.opposite().cell() uvw in the f.cell() uvw.
// the objective is to grow region where GPTransitionFunction makes each local frame compatible with the seed frame
void GPTransitionFunction::align_rot(Volume::Facet f, mat3x3 M_f_cell, mat3x3 M_opp_cell) {
    int best_mid = 0;
    double min_error = 1e20;
    for (ap.mid = 0; ap.mid < 24; ap.mid++) {
        mat3x3 M = ap.get_mat() * M_opp_cell - M_f_cell;
        double error = 0;//not M.sumsqr(): we only consider the restriction to f
        FOR(lv, 2) {
            vec3 e = f.vertex(lv + 1).pos() - f.vertex(0).pos();
            error += (M * e).norm2();
        }
        if (error < 1e-18) return;// fast-track (motly for identity)

        if (error < min_error) {
            best_mid = ap.mid;
            min_error = error;
        }
    }
    if (std::abs(min_error) > 1e-18) plop(min_error);

    ap.mid = best_mid;
}

GPTransitionFunction::GPTransitionFunction(Volume::Facet f, CellCornerAttribute<vec3>& U) {

    um_assert(f.opposite().active());

    Volume::Halfedge h[3] = { f.halfedge(0),f.halfedge(1),f.halfedge(2) };
    //auto opp = h[0].opposite_c();

    Volume::Halfedge opp = f.opposite().halfedge(0);
    int h_to = h[0].to();
    while ((opp.from() != h_to)) opp = opp.next();



    Volume::Halfedge opp_h[3] = { opp,opp.prev(),opp.next() };


    int best_mid = 0;
    double min_error = 1e20;
    for (ap.mid = 0; ap.mid < 24; ap.mid++) {
        mat3x3 M = ap.get_mat();
        double error = 0;//not M.sumsqr(): we only consider the restriction to f
        FOR(lh, 2) {
            vec3 e_c = U[h[(lh + 1) % 3].from_corner()] - U[h[lh].from_corner()];
            vec3 e_opp = U[opp_h[lh].from_corner()] - U[opp_h[lh].to_corner()];
            error += (e_c - M * e_opp).norm2();
        }
        if (error < min_error) {
            best_mid = ap.mid;
            min_error = error;
            if (error < 1e-18) {
                ap.mid = best_mid;
                t = U[h[0].from_corner()] - ap.get_mat() * U[opp_h[0].to_corner()];
                return;// fast-track (motly for identity)
            }
        }
    }
    plop(min_error);
    //um_assert(!"GPTransitionFunction not found on facet");
}

GPTransitionFunction GPTransitionFunction::inverted() {
    GPTransitionFunction res = *this;
    res.ap = ap.inverse();
    res.t = -(res.ap.get_mat() * t);
    return res;
}

GPTransitionFunction GPTransitionFunction::apply(GPTransitionFunction gp) {
    GPTransitionFunction res;
    res.ap = ap * gp.ap;
    res.t = ap.get_mat() * gp.t + t;
    return res;
}

void GPTransitionFunction::show() {
    std::cerr << "Rotation id : " << ap.mid << "\t  translation : " << t << std::endl;
}

vec3 GPTransitionFunction::apply(vec3 x) { return ap.get_mat() * x + t; }









void test_GPTransitionFunction(Tetrahedra& m, CellCornerAttribute<vec3>& U) {
    Trace::step("test uvw to jacobian");
    for (auto c : m.iter_cells()) {
        mat3x3 J = uvw_to_jacobian(c, U);
        FOR(lv, 3) {
            vec3 e_pos = c.vertex(lv + 1).pos() - c.vertex(0).pos();
            vec3 uvw_pos = U[c.corner(lv + 1)] - U[c.corner(0)];
            um_assert((J * e_pos - uvw_pos).norm() < 1e-10);
        }
    }
    Trace::step("test that transition function across a facet makes opposite compatible with current cell");
    for (auto f : m.iter_facets())if (f.opposite().active()) {
        GPTransitionFunction tr(f, U);
        FOR(lh, 3) {
            vec3 diff = U[f.halfedge(lh).from_corner()] - tr.apply(U[f.halfedge(lh).opposite_c().to_corner()]);
            um_assert(diff.norm2() < 1e-20);
        }
        if (!tr.ap.is_identity() || tr.t.norm2() > 2e-20)
            tr.show();
    }
    Trace::step("check composition of transition functions");
    FOR(it, 10) {
        GPTransitionFunction tr[10];
        FOR(i, 10) tr[i].ap.mid = rand() % 24;
        FOR(i, 10) FOR(d, 3) tr[i].t[d] = rand_range(-10., 10.);


        GPTransitionFunction compo = tr[0];
        FOR(i, 9) compo = tr[i + 1].apply(compo);

        vec3 test; FOR(d, 3) test[d] = rand_range(-10., 10.);
        vec3 res = test;
        for (int i = 9; i >= 0; i--) res = tr[i].apply(res);
        std::cerr << res << " ==?  " << compo.apply(test) << "  norm2 = " << (res - compo.apply(test)).norm2() << std::endl;
    }
    Trace::step("check inversion");
    FOR(it, 10) {
        GPTransitionFunction tr;
        tr.ap.mid = rand() % 24;
        FOR(d, 3) tr.t[d] = rand_range(-10., 10.);
        GPTransitionFunction compo = tr.apply(tr.inverted());
        vec3 test;
        FOR(d, 3) test[d] = rand_range(-10., 10.);
        std::cerr << test << " ==?  " << compo.apply(test) << "  norm2 = " << (test - compo.apply(test)).norm2() << std::endl;
    }
}
