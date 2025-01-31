#include <algo/volume/hex_smooth.h>

#include <ultimaille/all.h>
#include <framework/trace.h>

#include "toolbox.h"
#include "drop_attribute.h"

#include <iostream>


#include <algo/volume/hex_select.h>




void HexSmoother::show_scaled_jacobien(std::string name, double max_value) {
	m.connect();
	HexPointSelect badSJ(m, false);
	CellAttribute<double> SJ(m);
	ToolBox(m).eval_quality(SJ);
	double minSJ = 1e20;
	for (auto c : m.iter_cells()) minSJ = std::min(minSJ, SJ[c]);
	for (auto c : m.iter_cells()) if (SJ[c] < max_value) badSJ.set_cell(c, true);
	if (minSJ < max_value) {
		if (max_value < 1) badSJ.drop_cells(name + "badSJ4_", 4);
		badSJ.drop_cells(name + ((max_value < 1) ? "badSJ8_" : "SJ"), 8);
	}
}


using namespace UM::Linear;
void HexSmoother::smooth_LS(double constraint_weight) {
 	LeastSquares ls(3 * m.nverts());
	for (auto v : m.iter_vertices()) if (lock[v]) FOR(d, 3) ls.fix(3 * v + d, v.pos()[d]);
	for (auto h : m.iter_halfedges()) FOR(d, 3)  ls.add_to_energy(X(3 * h.from() + d) - X(3 * h.to() + d));

	constraint_weight = 100.;
	for (auto constraint : tan_constraints) {
		LinExpr eq = -constraint.n * constraint.P;
		FOR(d, 3) eq = eq + (constraint.n[d]*X(d + 3 * constraint.v));
		ls.add_to_energy(constraint_weight * eq);


		//FOR(d, 3) ls.add_to_energy(constraint_weight * (X(d + 3 * constraint.v)-constraint.P[d]));
		//ls.add_to_energy(constraint_weight * (constraint.n[0]*X(0 + 3 * constraint.v)+constraint.n[1]*X(1 + 3 * constraint.v)+constraint.n[2]*X(2 + 3 * constraint.v)-constraint.n * constraint.P));
	}
	ls.solve();
	for (auto v : m.iter_vertices()) FOR(d, 3) v.pos()[d] = ls.value(3 * v + d);
	DropVolume(m).apply("deformed vol");
}




void check_degenerated_topo(Hexahedra& hex) {
	EdgeGraph eg(hex);
	EdgeAttribute<bool> done(eg, false);

	PolyLine pl;
	EdgeAttribute<int> t(pl);
	pl.points = hex.points;
	for (auto h : hex.iter_halfedges()) {
		if (done[eg.edge_from_halfedge(h)])	continue;
		for (auto cir : h.iter_CCW_around_edge()) {
			done[eg.edge_from_halfedge(h)] = true;
			if (cir != h && cir.next().to() == h.next().to())
				FOR(i,2) {
					int s = pl.create_edges(1);
					if (i == 0) {
						t[s] = 0;
						pl.vert(s, 0) = h.from();
						pl.vert(s, 1) = h.to();
					} else {
						t[s] = 1;
						pl.vert(s, 0) = cir.from();
						pl.vert(s, 1) = cir.to();
					}
				}
		}
	}
	if (pl.nedges() > 0)
		Drop(pl, t)._force_radius(.1).apply_arrow("pl");
}

using Clock = std::chrono::high_resolution_clock;
using namespace std::literals::chrono_literals;
using namespace UM;

//double tet_volume(const vec3& A, const vec3& B, const vec3& C, const vec3& D) {
//	return ((B - A) * cross(C - A, D - A)) / 6.;
//}

inline double chi(double eps, double det) {
	if (det > 0)
		return (det + std::sqrt(eps * eps + det * det)) * .5;
	return .5 * eps * eps / (std::sqrt(eps * eps + det * det) - det);
}

inline double chi_deriv(double eps, double det) {
	return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}




void HexSmoother::smooth_elliptic(double constraint_weight) {

	um_assert(m.connected());
	check_degenerated_topo(m);
	double total_volume = 0;
	for (auto c : m.iter_cells()) {
		double v = Hexahedron(c).volume();
		total_volume += v;
	}
	um_assert(total_volume > 0);

	constexpr double theta = 1. / 1000.; // the energy is (1-theta)*(shape energy) + theta*(area energy)
	constexpr int bfgs_maxiter = 300; // max number of inner iterations
	constexpr int outer_maxiter = 10;  // max number of outer iterations
	constexpr double bfgs_threshold = 1e-5;
	constexpr double outer_threshold = 1e-3;



	CellCornerAttribute<mat<4, 3>> reference(m);    // desired tet geometry
	CellCornerAttribute<double>   volume(m);       // 3D triangle area
	std::vector<double> X(m.nverts() * 3, 0.); // optimization variables
	CellCornerAttribute<double>   detc(m);       // 3D triangle area

	FOR(v, m.nverts()) FOR(d, 3)  X[3 * v + d] = m.points[v][d];


	// precompute stuff for speedup
	std::vector<int> local_resp_h1;
	std::vector<int> local_resp_h2;
	std::vector<int> local_resp_h3;
	for (auto h1 : m.iter_halfedges()) {
		auto h2 = h1.opposite_f().next();
		auto h3 = h2.opposite_f().next();
		if (h1 < h2 || h1 < h3) continue;
		local_resp_h1.push_back(h1);
		local_resp_h2.push_back(h2);
		local_resp_h3.push_back(h3);
		if (h1 > 23) break;
	}

	std::vector<int> active_cells;
	for (int cell = 0; cell < m.ncells(); cell++) {
		bool has_free_vertex = false;
		FOR(lv, 8) has_free_vertex = has_free_vertex || !lock[m.vert(cell, lv)];
		if (has_free_vertex) active_cells.push_back(cell);
	}




	for (int cid = 0; cid < active_cells.size(); cid++) FOR(lh, 8) {
		int cell = active_cells[cid];
		auto h1 = Volume::Halfedge(m,24 * cell + local_resp_h1[lh]);
		auto h2 = Volume::Halfedge(m, 24 * cell + local_resp_h2[lh]);
		auto h3 = Volume::Halfedge(m, 24 * cell + local_resp_h3[lh]);


		double l = std::cbrt(total_volume / m.ncells()); // target hex edge length
		vec3 A = { 0,0,0 };
		vec3 B = { l,0,0 };
		vec3 C = { 0,l,0 };
		vec3 D = { 0,0,l };

		mat<3, 3> ST = { {B - A, C - A, D - A} };
		int c = h1.from_corner();
		reference[c] = mat<4, 3>{ { {-1,-1,-1},{1,0,0},{0,1,0},{0,0,1} } }*ST.invert_transpose();
		volume[c] = Tetrahedron(A, B, C, D).volume();
	}

	const auto getJ = [&](const std::vector<double>& X, Volume::Halfedge& h1, const int* verts)->mat<3, 3> { // get Jacobian matrix for tetrahedron t
		int c = h1.from_corner();
		mat<3, 3> J = {};
		for (int i : {0, 1, 2, 3})
			for (int d : {0, 1, 2})
				J[d] += reference[c][i] * X[verts[i] * 3 + d];
		return J;
		};

	auto starting_time = Clock::now();
	std::vector<SpinLock> spin_locks(X.size());

	double mindet = 0.;
	int ninverted = 0;
	for (int cid = 0; cid < active_cells.size(); cid++) FOR(lh, 8) {
		int cell = active_cells[cid];

		auto h1 = Volume::Halfedge(m, 24 * cell + local_resp_h1[lh]);
		auto h2 = Volume::Halfedge(m, 24 * cell + local_resp_h2[lh]);
		auto h3 = Volume::Halfedge(m, 24 * cell + local_resp_h3[lh]);

		const int verts[4] = { h1.from(),h1.to(),h2.to(),h3.to() };

		const mat<3, 3> J = getJ(X, h1, verts);
		double det = J.det();
		mindet = std::min(mindet, det);
		ninverted += (det <= 0);
	}

	constexpr double e0 = 1e-4;
	double eps = mindet > 0 ? e0 : std::sqrt(e0 * e0 + 0.004 * mindet * mindet);
	std::cerr << "Untangling nverts: " << m.nverts() << std::endl;

	for (int iter = 0; iter < outer_maxiter; iter++) {
		std::cerr << "Outer iteration #" << iter << ", inverted corners: " << ninverted << std::endl;

		const STLBFGS::func_grad_eval func = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
			std::fill(G.begin(), G.end(), 0);
			F = 0;
			ninverted = 0;
			mindet = std::numeric_limits<double>::max();
			double local_min_det = mindet;
			int local_ninverted = ninverted;
			double local_F = F;



			int ONE = std::numeric_limits<int>::max();
#pragma omp parallel for reduction(+:local_F) 
			for (int constraint_id = 0; constraint_id < tan_constraints.size(); constraint_id++) {
				auto constraint = tan_constraints[constraint_id];
				LinExpr eq = constraint.n[0] * Linear::X(3 * constraint.v+0) + constraint.n[1] * Linear::X(3 * constraint.v+1) + constraint.n[2] * Linear::X(3 * constraint.v+2) + 
					
					constraint.n * -constraint.P * Linear::X(ONE) ;
				eq = constraint_weight * eq;
				FOR(i, eq.size()) FOR(j, eq.size()) {
					// participation to the function 
					double c_i = eq[i].value * (eq[i].index == ONE ? 1. : X[eq[i].index]);
					double c_j = eq[j].value * (eq[j].index == ONE ? 1. : X[eq[j].index]);
					local_F += c_i * c_j;
					// participation to the gradiant
					if (eq[i].index == ONE || lock[eq[i].index / 3]) continue;
					spin_locks[eq[i].index].lock();
					G[eq[i].index] += eq[i].value * eq[j].value * (eq[j].index == ONE ? 1. : X[eq[j].index]);
					spin_locks[eq[i].index].unlock();
				}

			}



#pragma omp parallel for reduction(+:local_F) reduction(min:local_min_det) reduction(+:local_ninverted)
			for (int cid = 0; cid < active_cells.size(); cid++) FOR(lh, 8) {
				int cell = active_cells[cid];

				Volume::Halfedge h1(m, 24 * cell + local_resp_h1[lh]);
				Volume::Halfedge h2(m, 24 * cell + local_resp_h2[lh]);
				Volume::Halfedge h3(m, 24 * cell + local_resp_h3[lh]);


				const int verts[4] = { h1.from(),h1.to(),h2.to(),h3.to() };

				const mat<3, 3> J = getJ(X, h1, verts);
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
				detc[h1.from_corner()] = det;
				local_min_det = std::min(local_min_det, det);
				local_ninverted += (det <= 0);

				double c1 = chi(eps, det);
				double c3 = chi_deriv(eps, det);
				double c2 = pow(c1, 2. / 3.);

				double f = J.sumsqr() / (3. * c2);
				double g = (1. + det * det) / (2. * c1);


				int corner_1 = h1.from_corner();

				local_F += ((1 - theta) * f + theta * g) * volume[corner_1];


				for (int d : {0, 1, 2}) {
					const vec3& a = J[d];
					const vec3& b = K[d];

					const vec3 dfda = (a * 2.) / (3. * c2) - b * ((2. * f * c3) / (3. * c1));
					const vec3 dgda = b * ((det - g * c3) / c1);

					for (int i : {0, 1, 2, 3}) {
						const int v = verts[i];
						if (lock[v]) continue;
						spin_locks[v * 3 + d].lock();
						G[v * 3 + d] += (dfda * (1. - theta) + dgda * theta) * reference[corner_1][i] * volume[corner_1];
						spin_locks[v * 3 + d].unlock();
					}
				}
			}
			mindet = local_min_det;
			ninverted = local_ninverted;
			F = local_F;
			};


		double E_prev, E;
		std::vector<double> trash(X.size());
		func(X, E_prev, trash);

		STLBFGS::Optimizer opt{ func };
		opt.ftol = opt.gtol = bfgs_threshold;
		opt.maxiter = bfgs_maxiter;
		opt.run(X);

		func(X, E, trash);
		std::cerr << "E: " << E_prev << " --> " << E << ", eps: " << eps << ", min det: " << mindet << std::endl;

		const double sigma = std::max(1. - E / E_prev, 1e-1);
		double mu = (1 - sigma) * chi(eps, mindet);
		if (mindet < mu)
			eps = 2 * std::sqrt(mu * (mu - mindet));
		else eps = 1e-10;
		if (mindet > 0 && std::abs(E_prev - E) / E < outer_threshold) break;
	}

	std::cerr << "Running time: " << (Clock::now() - starting_time) / 1.s << " seconds" << std::endl;
	for (int v : m.iter_vertices())
		m.points[v] = { X[3 * v + 0], X[3 * v + 1], X[3 * v + 2] };
}


