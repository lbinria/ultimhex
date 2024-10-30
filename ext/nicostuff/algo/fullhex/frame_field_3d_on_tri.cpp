
#include <drop_attribute.h>
#include <framework/trace.h>
#include <fullhex/gp_basic.h>
#include <fullhex/frame_field_3d_on_tri.h>

#include <volume/frame3D.h>


using namespace UM::Linear;
FF3D_on_tri::FF3D_on_tri(Triangles& m) :m(m), r(m), J(m) {}

void FF3D_on_tri::show_surface_rot(std::string name, double scale, bool with_flag) {
    double ave = ToolBox(m).ave_edge_size();
    Hexahedra outm;
    outm.create_cells(m.nfacets());
    outm.points.create_points(8 * m.nfacets());
    CellFacetAttribute<int> flag(outm, 0);
    FOR(f, outm.nfacets())flag[f] = f % 6;
    for (auto f : m.iter_facets()) {
        FOR(v, 8) outm.vert(f, v) = 8 * f + v;
        FOR(di, 2)FOR(dj, 2)FOR(dk, 2) outm.points[8 * f + di + 2 * dj + 4 * dk] = Triangle3(f).bary_verts() + scale * .15 * ave
            * J[f].transpose() /*J is a rot => transpose == fast_invert*/ * vec3(2 * di - 1, 2 * dj - 1, 2 * dk - 1);
    }
    if (with_flag)	Drop(outm, flag).apply(name + "_cubes");
    else			DropVolume(outm).apply(name + "_cubes");
}

void FF3D_on_tri::init_rot_by_FF(double normal_coeff) {
    //DropSurface(m).apply("in surface");
    double smooth_coeff = 1;
    LeastSquares ls(11 * m.nfacets(), true);

    Trace::step("Construct matrix");
    for (auto h : m.iter_halfedges()) {
        auto opp = h.opposite();
        if (!opp.active()) continue;
        FOR(i, 9) ls.add_to_energy(smooth_coeff * (X(11 * opp.facet() + i) - X(11 * h.facet() + i)));
    }


    // boundary condition enforced by barrier equations
    for (auto f : m.iter_facets()) {
        J[f] = mat3x3::identity();
        vec3 n = Triangle3(f).normal();
        vec3 cr = -cross(n, vec3(0, 0, 1));
        if (cr.norm2() > 1e-10) {
            double angle = std::atan2(cr.norm(), n[2]);
            Quaternion q = vec3_to_quat(angle * cr.normalized());
            J[f] = q.rotation_matrix();
        }

        SH4 sh0, sh4, sh8;
        sh4[4] = std::sqrt(7. / 12.);
        sh0[0] = std::sqrt(5. / 12.);
        sh8[8] = std::sqrt(5. / 12.);
        vec3 xyz = mat_to_euler(normalize_columns(J[f]));
        sh4.euler_rot(xyz);
        sh0.euler_rot(xyz);
        sh8.euler_rot(xyz);

        FOR(i, 9) ls.add_to_energy(normal_coeff * (
            X(11 * f + i)
            + sh0[i] * X(11 * f + 9)
            + sh8[i] * X(11 * f + 10)
            - sh4[i]
            ));
    }

    Trace::step("Solve SH");
    ls.solve();



    Trace::step("Project SH to B");
    // convert spherical harmonic coefficients to a rotation
    {
#pragma omp parallel for
        FOR(f, m.nfacets()) {
            SH4 fv;
            FOR(i, 9) fv[i] = ls.value(f * 11 + i);
            J[f] = fv.project_mat(1e-3, 1e-5, nullptr);//fv.project_mat(1e-3, 1e-5, nullptr);
            J[f] = J[f].transpose();
        }
    }
}


void FF3D_on_tri::init_rotation_field_polycube() {
    for (auto f : m.iter_facets()) {
        vec3 n = Triangle3(f).normal();
        int d = 0;
        FOR(i, 2) if (std::abs(n[i + 1]) > std::abs(n[d])) d = i + 1;
        vec3 obj = mat3x3::identity()[d];
        if (n[d] < 0)obj *= -1;
        r[f] = cross(n, obj);
        if (r[f].norm2() < 1e-10) {
            r[f] = vec3(0, 0, 0);
            continue;
        }
        double angle = std::atan2(r[f].norm(), n * obj);
        r[f] = angle * r[f].normalized();
    }
    // output J
    FOR(f, m.nfacets()) J[f] = vec3_to_quat(r[f]).rotation_matrix();

}

void FF3D_on_tri::rigid_rotation_mono_25D() {
    // find best global rotation
    std::vector<vec3> normal(m.nfacets());
    std::vector<double> w(m.nfacets());
    for (auto f : m.iter_facets()) {
        Triangle3 tri = Triangle3(f);
        normal[f] = tri.normal();
        w[f] = tri.unsigned_area();
    }
    mat3x3 R = Frame::representative_frame(normal, w);
    // rotate the object
    for (auto v : m.iter_vertices())
        v.pos() = R * v.pos();
    plop(R);


    // permute axis to set z as stable coordinate
    double z_score[3] = { 0,0,0 };
    for (auto f : m.iter_facets()) {
        vec3 n = Triangle3(f).normal();
        FOR(d, 3) z_score[d] += std::cos(n * R[d]*2.*M_PI);
    }

    if (z_score[0] > z_score[2] && z_score[0] > z_score[1])
        for (auto v : m.iter_vertices())
            v.pos() = vec3(v.pos()[1], v.pos()[2], v.pos()[0]);
    if (z_score[1] > z_score[2] && z_score[1] > z_score[0])
        for (auto v : m.iter_vertices())
            v.pos() = vec3(v.pos()[2], v.pos()[0], v.pos()[1]);
}


void FF3D_on_tri::init_rotation_field_mono_25D(bool smooth_transition) {

    // addapt rotation to normal
    for (auto f : m.iter_facets()) {
        vec3 n = Triangle3(f).normal();
        r[f] = cross(n, vec3(0, 0, 1));
        double angle = std::atan2(r[f].norm(), n[2]);
        if (std::abs(n[2]) > .99) {
            r[f] = vec3(0, 0, 0);
            continue;
        }

        if (smooth_transition)
            angle = std::sin(angle * 4.) / 4.;
        else {
            if (angle < M_PI / 4.);
            else if (angle < 3. * M_PI / 4.) angle = angle - M_PI / 2.;
            else angle = angle - M_PI;
        }

        r[f] = angle * r[f].normalized();
    }
    // output J
    FOR(f, m.nfacets()) J[f] = vec3_to_quat(r[f]).rotation_matrix();

}



void FF3D_on_tri::smooth_rotation_field(double fit_coeff) {
    LeastSquares ls(m.nfacets() * 3);

    // data fitting
    for (auto f : m.iter_facets()) {
        if (Triangle3(f).unsigned_area() < 1e-10) continue;
        double angle = r[f].norm();

        vec3 B[3];
        B[0] = Triangle3(f).normal();
        B[1] = r[f];
        if (angle < 1e-10) {
            B[1] = cross(vec3(1, 0, 0), B[0]);
            if (B[1].norm2() < .1)
                B[1] = cross(vec3(0, 1, 0), B[0]);
        }
        B[1].normalize();
        B[2] = cross(B[0], B[1]);

        FOR(b, 3) {
            double rhs = 0;
            if (b == 1 && angle > 1e-10) rhs = angle;
            double w = 1;
            if (b == 0) w = .01;
            ls.add_to_energy(w * fit_coeff * (
                B[b][0] * X(3 * f + 0) + B[b][1] * X(3 * f + 1) + B[b][2] * X(3 * f + 2) - rhs)
            );
        }
    }
    // smooth term
    for (auto h : m.iter_halfedges()) {
        um_assert(!h.on_boundary());
        auto f = h.facet();
        auto fo = h.opposite().facet();
        FOR(d, 3) ls.add_to_energy(X(3 * f + d) - X(3 * fo + d));
    }

    ls.solve();
    for (auto f : m.iter_facets()) FOR(d, 3) r[f][d] = ls.X[3 * f + d];

    // output J
    FOR(f, m.nfacets()) J[f] = vec3_to_quat(r[f]).rotation_matrix();

}
