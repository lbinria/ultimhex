#include <framework/mesh_API.h>
#include <fullhex/gp_basic.h>
#include <fullhex/frame_field_3d_on_tri.h>
using namespace UM::Linear;
const std::vector<std::string> FF3D_cmd::algorithms = {
    "spherical_harmonics_tet","spherical_harmonics_tri","polycube_tri","25d_tri"
};

void FF3D_cmd::apply(std::string algo, double constraint_weight) {
    Tetrahedra& tet = mesh3d.m;
    CellAttribute<mat3x3> J("J",mesh3d.attr,mesh3d.m);
    FF3D ff(tet, J);

    CellFacetAttribute<bool> lock_normal(tet, false);

    if (constraint_weight < .0001)
        Trace::abort("Without constraints, the frame field would just be constant... not implemented yet");


    if (algo == "spherical_harmonics_tet") {
        for (auto f : tet.iter_facets()) lock_normal[f] = f.on_boundary();
        ff.apply_ls(lock_normal, 100. * constraint_weight);
    }

    else {
        TetBoundary bound(tet);
        FF3D_on_tri ff_tri(bound.tri);

        if (algo == "spherical_harmonics_tri") {
            ff_tri.init_rot_by_FF(100. * constraint_weight);
            CellAttribute<bool> lock_cell(tet, false);
            for (auto f : bound.tri.iter_facets()) {
                J[bound.tet_facet(f).cell()] = ff_tri.J[f];
                lock_cell[bound.tet_facet(f).cell()] = true;
            }
            ff.apply_ls(lock_normal, 100, &lock_cell);
        }
        else {
            if (algo == "polycube_tri") {
                ff_tri.init_rotation_field_polycube();
                ff_tri.smooth_rotation_field(100. * constraint_weight);
            }
            if (algo == "25d_tri") {
                ff_tri.init_rotation_field_mono_25D(false);
                ff_tri.smooth_rotation_field(100. * constraint_weight);
            }
            // propagate inside without singularities (interpolates r instead of J)
            CellAttribute<vec3> r(tet);
            FOR(dim, 3) {
                LeastSquares ls(tet.ncells());
                for (auto f : bound.tri.iter_facets())
                    ls.fix(bound.tet_facet(f).cell(), ff_tri.r[f][dim]);
                for (auto f : tet.iter_facets()) if (!f.on_boundary())
                    ls.add_to_energy(X(f.cell()) - X(f.opposite().cell()));
                ls.solve();
                for (auto c : tet.iter_cells())
                    r[c][dim] = ls.value(c);
            }
            for (auto c : tet.iter_cells()) J[c] = vec3_to_quat(r[c]).rotation_matrix();
        }
    }


    ff.show_cubes("cube_" + algo);
    ff.show_singularity_graph("singu_" + algo);
    CellCornerAttribute<vec3> ff_export(tet);

    //mesh3d.add(J, "J");
}