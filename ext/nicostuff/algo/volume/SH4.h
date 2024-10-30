
#pragma once

#include <ultimaille/all.h>
#include <basic.h>

using namespace UM;

struct SH4 {

	vec<9> coeff;

	SH4() {
		FOR(i, 9)  coeff[i] = 0.;
	}

	SH4(const vec<9>& rhs) : coeff(rhs) {
	}

	SH4(double* fv) {
		FOR(i, 9)  coeff[i] = fv[i];
	}

	SH4(
		double x0, double x1, double x2,
		double x3, double x4, double x5,
		double x6, double x7, double x8
	) {
		coeff[0] = x0; coeff[1] = x1; coeff[2] = x2;
		coeff[3] = x3; coeff[4] = x4; coeff[5] = x5;
		coeff[6] = x6; coeff[7] = x7; coeff[8] = x8;
	}

	double& operator[](int i) {
		um_assert(i < 9);
		return coeff[i];
	}

	double norm() const {
		return coeff.norm();
	}

	double operator *(const SH4& other) const {
		return coeff * other.coeff;
	}

	SH4 operator -(const SH4& other) const {
		return SH4(coeff - other.coeff);
	}

	SH4 operator *(double s) const {
		return SH4(s * coeff);
	}

	SH4 operator /(double s) const {
		return SH4(coeff / s);
	}

	SH4 operator +(const SH4& v) const {
		return SH4(coeff + v.coeff);
	}

	static double basis(int id, const vec3& v);

	double value(const vec3& v) const {
		double res = 0;
		FOR(i, 9)res += coeff[i] * basis(i, v);
		return res;
	}


	void Rz(double alpha);
	void Ry(double alpha);
	void Rx(double alpha);

	void euler_rot(const vec3& rot_vec) {
		Rx(rot_vec[0]);
		Ry(rot_vec[1]);
		Rz(rot_vec[2]);
	}

	SH4 Ex() const;
	SH4 Ey() const;
	SH4 Ez() const;

	static SH4 rest_frame() {
		return SH4(0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, std::sqrt(5. / 12.));
	}

	mat<3, 3> project_mat(double grad_threshold = 1e-3, double dot_threshold = 1e-5, vec3* euler_prev = nullptr);

};
