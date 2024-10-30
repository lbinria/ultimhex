#ifndef FF_ON_TRI__H__
#define FF_ON_TRI__H__


#include <framework/trace.h>
#include <fullhex/gp_basic.h>


struct FF3D_on_tri {
    FF3D_on_tri(Triangles& m);

    void show_surface_rot(std::string name = "rot", double scale = 1., bool with_flag = true);

    void init_rot_by_FF(double normal_coeff = 100.);


    void init_rotation_field_polycube();
    
    
    void rigid_rotation_mono_25D();
    void init_rotation_field_mono_25D(bool smooth_transition = true);
    void smooth_rotation_field(double fit_coeff = .01);
    Triangles& m;
    FacetAttribute<vec3> r;
    FacetAttribute<mat3x3> J;
};
#endif