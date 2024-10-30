#ifndef DROP__GLYPH__H__
#define DROP__GLYPH__H__


#include <ultimaille/all.h>
#include <basic.h>

#include <algo/drop_mesh.h>
#include <algo/drop_attribute.h>




void drop_point(vec3 P, std::string name);
void drop_arrow(vec3 from, vec3 to, std::string name);
void drop_triangle(Triangle3 t, std::string name);


#endif