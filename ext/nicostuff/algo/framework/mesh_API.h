#ifndef MESH__API__H__
#define MESH__API__H__
#include <ultimaille/all.h>
#include <create_shape.h>
#include <framework/trace.h>
#include <framework/meta.h>
#include "toolbox.h"
#include "drop_attribute.h"

#include <iostream>
#include <fullhex/gp_basic.h>
#include <fullhex/frame_field_3d_on_tri.h>


inline std::string serialize_string_vector(const std::vector<std::string>& array) {
    std::string res;
    for (int i = 0;; i++) {
        res += array[i];
        if (i + 1 < array.size()) res += ",";
        else return res;
    }
}

struct FF3D_cmd {
    inline FF3D_cmd(Meta::Mesh3d<Tetrahedra>& mesh3d) : mesh3d(mesh3d) {}
    void apply(std::string algo, double constraint_weight);
    Meta::Mesh3d<Tetrahedra>& mesh3d;
    static const std::vector<std::string> algorithms;
};




#endif