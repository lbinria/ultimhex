#include <fullhex/untangle_tetrahedra.h>
#include <ultimaille/all.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <chrono>
#include <basic.h>
#include <toolbox_tetrahedra.h>
#include <toolbox_hexahedra.h>

    const vec3 unit_regular_tet[4] = {
            vec3{ 0,   0, 0},vec3{ -1,   0, 0},vec3{ -.5, -.5,  1. / std::sqrt(2.) },vec3{ -.5,  .5,  1. / std::sqrt(2.) }
        };


using Clock = std::chrono::high_resolution_clock;
using namespace std::literals::chrono_literals;
using namespace UM;
using namespace UM::Linear;

inline double chi(double eps, double det) {
    if (det > 0) return (det + std::sqrt(eps * eps + det * det)) * .5;
    return .5 * eps * eps / (std::sqrt(eps * eps + det * det) - det);
}

inline double chi_deriv(double eps, double det) {
    return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}

double quality(const vec3& A, const vec3& B, const vec3& C, const vec3& D) {
    //constexpr
    mat3x3 Sinv = { {{1, -1. / sqrt(3.), -1. / sqrt(6.)}, {0, 2 / sqrt(3.), -1. / sqrt(6.)}, {0, 0, sqrt(1.5)}} };
    const mat3x3 J = mat3x3{ {B - A,C - A,D - A} }.transpose() * Sinv;
    const double det = J.det();
    um_assert(det > 0);
    const mat3x3 Ji = J.invert();
    auto [U1, D1, V1] = svd3x3(J);
    auto [U2, D2, V2] = svd3x3(Ji);
    return D1[0][0] * D2[0][0];  // L2 norm of S1 times L2 norm of S2; it is possible to use Frobenius condition number as J.norm()*Ji.norm();
}


// low-level API for untangling tetrahedra
// 
constexpr double reference_quality_threshold = 7.; // replace all reference tetrahedra with condition number exceeding the threshold with equilateral tet reference
void untangle_tetrahedra(Tetrahedra& m, CellAttribute<mat<4, 3>>& reference, CellAttribute<double>& volume, std::vector<LinExpr>* ls_matrix_ptr = NULL, std::vector<LinExpr>* hard_constraints_matrix_ptr = NULL) {
    //constexpr double theta = 1e-3; // the energy is (1-theta)*(shape energy) + theta*(volume energy)
    constexpr double theta = 1e-2; // the energy is (1-theta)*(shape energy) + theta*(volume energy)
    constexpr int bfgs_maxiter = 3000; // max number of inner iterations
    constexpr int outer_maxiter = 150;  // max number of outer iterations
    constexpr double bfgs_threshold = 1e-8;
    constexpr double outer_threshold = 1e-3;
    
    int nX = 3 * m.nverts();
    NullSpaceBuilder rb(nX + 1, true);
    std::vector<LinExpr> empty_constraints;
    std::vector<LinExpr>& hard_constraints_matrix = (hard_constraints_matrix_ptr == NULL) ? empty_constraints : *hard_constraints_matrix_ptr;

    FOR(line_id, hard_constraints_matrix.size()) {
        SparseVector eq = hard_constraints_matrix[line_id];
        for (auto& c : eq) if (c.index == -1) c.index = nX;
        rb.add_constraint(eq);
    }
    CRSMatrix M, Mt;
    M = rb.to_crs();
    Mt = M.transpose();

    int nY = Mt.nrows();
    std::vector<double> Y(nY, 0.);
    Y[nY-1] = 1.;

    
    // init Y from vertices pos 
    for (int row = 0; row < nX; ++row) {
            auto it = M.iter_row(row);
            const SparseElement &e = *it.begin();
            if (it.end()-it.begin() == 1 && e.index < nY-1) {
                um_assert(std::abs(e.value-1.)<1e-10);
                Y[e.index] = m.points[row/3][row%3];
            }
        }


    const auto getJ = [&reference, &m](const std::vector<double>& X, int t)->mat<3, 3> { // get Jacobian matrix for tetrahedron t
        mat<3, 3> J = {};
        for (int i : {0, 1, 2, 3}) FOR(d, 3) J[d] += reference[t][i] * X[m.vert(t, i) * 3 + d];
        return J;
        };

    auto starting_time = Clock::now();
    std::vector<SpinLock> spin_locks(nX);

    double mindet = 0.;
    int ninverted = 0;
    {
        // init mindet and ninverted
        std::vector<double> X = M * Y;
        for (auto t : m.iter_cells()) {
            const mat<3, 3> J = getJ(X, t);
            double det = J.det();
            mindet = std::min(mindet, det);
            ninverted += (det <= 0);
        }
    }

    constexpr double e0 = 1e-4;
    double eps = mindet > 0 ? e0 : std::sqrt(e0 * e0 + 0.004 * mindet * mindet); // heuristic choice of the first value for the continuation parameter

    FOR(iter, outer_maxiter) {
        std::cerr << "Outer iteration #" << iter << ", inverted tetrahedra: " << ninverted << std::endl;

        const STLBFGS::func_grad_eval func = [&](const std::vector<double>& y, double& F, std::vector<double>& gy) {
            std::vector<double> x = M * y;
            std::vector<double> gx(nX+1, 0);
            //std::fill(G.begin(), G.end(), 0);
            int omp_ninverted = 0;
            double omp_mindet = std::numeric_limits<double>::max();
            double omp_F = 0;

            if (ls_matrix_ptr != NULL) {
                std::vector<LinExpr>& ls_matrix = *ls_matrix_ptr;
                //plop(ls_matrix .size());
                
#pragma omp parallel for reduction(+:omp_F) 
                FOR(line_id, ls_matrix.size()) {
                    //plop(line_id);
                    LinExpr eq = ls_matrix[line_id];
                    //plop(eq);
                    FOR(i, eq.size()) FOR(j, eq.size()) {
                        //plop(j);
                        // participation to the function 
                        double c_i = eq[i].value * (eq[i].index == -1 ? 1. : x[eq[i].index]);
                        double c_j = eq[j].value * (eq[j].index == -1 ? 1. : x[eq[j].index]);
                        omp_F += c_i * c_j;
                        // participation to the gradiant
                        //if (eq[i].index == -1 || lock[eq[i].index / 3]) continue;
                        if (eq[i].index!=-1){
                            spin_locks[eq[i].index].lock();
                            gx[eq[i].index] += eq[i].value * eq[j].value * (eq[j].index == -1 ? 1. : x[eq[j].index]);
                            spin_locks[eq[i].index].unlock();
                        }
                    }
                }
            }


#pragma omp parallel for reduction(+:omp_F) reduction(min:omp_mindet) reduction(+:omp_ninverted)
            FOR(t, m.ncells()) {
                const mat<3, 3> J = getJ(x, t);
                const mat<3, 3> K = { // dual basis
                    {{
                         J[1].y * J[2].z - J[1].z * J[2].y,
                         J[1].z * J[2].x - J[1].x * J[2].z,
                         J[1].x * J[2].y - J[1].y * J[2].x
                     },
                    {
                        J[0].z * J[2].y - J[0].y * J[2].z,
                        J[0].x * J[2].z - J[0].z * J[2].x,
                        J[0].y * J[2].x - J[0].x * J[2].y
                    },
                    {
                        J[0].y * J[1].z - J[0].z * J[1].y,
                        J[0].z * J[1].x - J[0].x * J[1].z,
                        J[0].x * J[1].y - J[0].y * J[1].x
                    }}
                };
                const double det = J[0] * K[0];
                omp_mindet = std::min(omp_mindet, det);
                omp_ninverted += (det <= 0);

                double c1 = chi(eps, det);
                double c3 = chi_deriv(eps, det);
                double c2 = pow(c1, 2. / 3.);

                double f = J.sumsqr() / (3. * c2);
                double g = (1. + det * det) / (2. * c1);
                omp_F += ((1 - theta) * f + theta * g) * volume[t];

                FOR(d,3) {
                    const vec3& a = J[d];
                    const vec3& b = K[d];

                    const vec3 dfda = (a * 2.) / (3. * c2) - b * ((2. * f * c3) / (3. * c1));
                    const vec3 dgda = b * ((det - g * c3) / c1);

                    FOR(i,4) {
                        const int v = m.vert(t, i);
                        //if (lock[v]) continue;
                        spin_locks[v * 3 + d].lock();
                        gx[v * 3 + d] += (dfda * (1. - theta) + dgda * theta) * reference[t][i] * volume[t];
                        spin_locks[v * 3 + d].unlock();
                    }
                }
            }
            ninverted = omp_ninverted;
            mindet = omp_mindet;
            F = omp_F;
            gy = Mt * gx;
            gy.back() = 0.;
            };


        double E_prev, E;
        std::vector<double> trash(Y.size());
        func(Y, E_prev, trash);

        STLBFGS::Optimizer opt(func, 10, true);
        opt.ftol = opt.gtol = bfgs_threshold;
        opt.maxiter = bfgs_maxiter;
        opt.run(Y);

        func(Y, E, trash);
        std::cerr << "f: " << E_prev << " --> " << E << ", eps: " << eps << ", min det: " << mindet << std::endl;

        const double sigma = std::max(1. - E / E_prev, 1e-1);
        double mu = (1 - sigma) * chi(eps, mindet);
        if (mindet < mu)
            eps = 2 * std::sqrt(mu * (mu - mindet));
        else eps = 1e-10;
        if (mindet > 0 && std::abs(E_prev - E) / E < outer_threshold) break;
    }

    std::cerr << "Running time: " << (Clock::now() - starting_time) / 1.s << " seconds" << std::endl;
    if (mindet > 0)     std::cerr << "SUCCESS" << std::endl;
    else                std::cerr << "FAIL" << std::endl;
    

    std::vector<double> X = M * Y;
    for (auto v : m.iter_vertices()) FOR(d, 3) v.pos()[d] = X[3 * v + d];
}



 
    UntanglerTet::UntanglerTet(Tetrahedra& m) :m(m), lock(m), reference(m), volume(m) {
        // lock boundary
        for (auto f : m.iter_facets()) if (f.on_boundary()) FOR(lv, 3) lock[f.vertex(lv)] = true;

        // init reference to regular tets
        double vol_ave = ToolBox(m).get_volume() / double(m.ncells());
        double a = std::cbrt(vol_ave * 6. * std::sqrt(2.));
        for (auto t : m.iter_cells()) {
            volume[t] = vol_ave;
            vec3 tet[4] ;            
            FOR(lv,4) tet[lv] = a * unit_regular_tet[lv]; 
            reference[t] = mat<4, 3>{ { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} } }*mat<3, 3>{{tet[1] - tet[0], tet[2] - tet[0], tet[3] - tet[0]}}.invert_transpose();
        }

        minq = +1e10;
        maxq = -1e10;

    }


    void UntanglerTet::check_obvious_overconstraints() {
        for (auto c : m.iter_cells()) {
            int nlock = 0;
            FOR(lv, 4) nlock += lock[c.vertex(lv)];
            um_assert(nlock < 4 || Tetrahedron(c).volume()>0);
        }
    }

    void UntanglerTet::set_tet_shape_objective(int tet_id,Tetrahedron tet) {
        volume[tet_id] = tet.volume();
        um_assert(volume[tet_id] > 0);
            double qual = quality(tet[0], tet[1], tet[2], tet[3]);
            minq = std::min(minq, qual);
            maxq = std::max(maxq, qual);
            if (qual > reference_quality_threshold)
            { // do we need to repair the reference?
                double a = std::cbrt(volume[tet_id] * 6. * std::sqrt(2.));
                FOR(lv, 4) tet[lv] = a * unit_regular_tet[lv];
            }
            reference[tet_id] = mat<4, 3>{ { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} } }*mat<3, 3>{{tet[1] - tet[0], tet[2] - tet[0], tet[3] - tet[0]}}.invert_transpose();
    }

    void UntanglerTet::set_shape_objective(Tetrahedra& ref) {
        for (auto t : ref.iter_cells()) 
            set_tet_shape_objective(t, Tetrahedron(t));
        
        std::cerr << "Input tetra condition number varies from " << minq << " to " << maxq << ", repair threshold was chosen at " << reference_quality_threshold << std::endl;
        //double minq = +1e10;
        //double maxq = -1e10;
        //for (auto t : ref.iter_cells()) {
        //    volume[t] = Tetrahedron(t).volume();
        //    um_assert(volume[t] > 0);
        //    vec3 tet[4];
        //    FOR(lv, 4) tet[lv] = t.vertex(lv).pos();
        //    double qual = quality(tet[0], tet[1], tet[2], tet[3]);
        //    minq = std::min(minq, qual);
        //    maxq = std::max(maxq, qual);
        //    if (qual > reference_quality_threshold)
        //    { // do we need to repair the reference?
        //        double a = std::cbrt(volume[t] * 6. * std::sqrt(2.));
        //    FOR(lv,4) tet[lv] = a * unit_regular_tet[lv]; 
        //    }
        //    reference[t] = mat<4, 3>{ { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} } }*mat<3, 3>{{tet[1] - tet[0], tet[2] - tet[0], tet[3] - tet[0]}}.invert_transpose();
        //}
        //std::cerr << "Input tetra condition number varies from " << minq << " to " << maxq << ", repair threshold was chosen at " << reference_quality_threshold << std::endl;
    }



    void UntanglerTet::normalize() {
        auto box = ToolBox(m.points).bbox();
        translation = box.center();
        scale = 10. / (box.max - box.min).norm();
        inv_scale = 1. / scale;
        if (ToolBox(m).get_volume() < 0) scale *= -1;
        for (auto c : m.iter_cells()) reference[c] = inv_scale * reference[c];
        for (auto c : m.iter_cells()) volume[c] *= scale * scale * scale;
        for (auto v : m.iter_vertices()) v.pos() -= translation;
        for (auto v : m.iter_vertices()) v.pos() *= scale;
    }

    void UntanglerTet::un_normalize() {
        for (auto v : m.iter_vertices()) v.pos() *= inv_scale;
        for (auto v : m.iter_vertices()) v.pos() += translation;
        for (auto c : m.iter_cells()) reference[c] = scale * reference[c];
        for (auto c : m.iter_cells()) volume[c] *= inv_scale * inv_scale * inv_scale;
    }

    void UntanglerTet::apply(bool need_normalize ) {
        check_obvious_overconstraints();
        if (need_normalize ) normalize();


        // test with hard constraints
        for (auto v : m.iter_vertices()) if (lock[v]) FOR(d, 3) hard_constraints_matrix.push_back(X(3 * v + d) - v.pos()[d]);
        
        
        // test with ls constraints
        //for (auto v : m.iter_vertices()) if (lock[v]) FOR(d, 3) ls_matrix.push_back(X(3 * v + d) - v.pos()[d]);

        untangle_tetrahedra(m,  reference, volume, ls_matrix.empty()? NULL:&ls_matrix, hard_constraints_matrix.empty()?NULL:&hard_constraints_matrix);
        if (need_normalize ) un_normalize();
    }

