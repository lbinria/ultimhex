#ifndef UNTANGLE_TET__H__
#define UNTANGLE_TET__H__

#include <ultimaille/all.h>
#include <basic.h>



using namespace UM;
using namespace UM::Linear;



/***
* UntanglerTet untangles a tet mesh
*     - supports extra energy terms as least squares energy
*     - supports linear constraints by DOF elimintation
**/
struct UntanglerTet {
    UntanglerTet(Tetrahedra& m) ;

    void check_obvious_overconstraints();

    void set_tet_shape_objective(int tet_id, Tetrahedron tet);
    void set_shape_objective(Tetrahedra& ref) ;



    void normalize() ;

    void un_normalize();

    void apply(bool need_normalize = true) ;


    std::vector<LinExpr>        ls_matrix;
    std::vector<LinExpr>        hard_constraints_matrix;
    CellAttribute<mat<4, 3>>    reference;      // desired tet geometry
    CellAttribute<double>       volume;         // could be derived from reference, maybe later
    PointAttribute<bool>        lock;           // locking some vertices
    vec3 translation;                           // move and scale the mesh to limit numerical issues
    double scale;
    double inv_scale;
    Tetrahedra& m;

    double minq ;
    double maxq ;

};

#endif