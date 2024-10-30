

#pragma once
#include <ultimaille/all.h>
#include <cmath>

#include "SH4.h" // to compute average_frame(vector<mat<3,3>>& data) and representative_frame

using namespace UM;


Quaternion align_with_uv(vec3 u, vec3 v);
Quaternion align_with_v(Quaternion const& q, vec3 v, int& axis);
vec3 quat_to_vec3(Quaternion const& q);
Quaternion vec3_to_quat(vec3 const& v);


const vec3 label_to_normal[6] = {
	vec3(1, 0, 0), vec3(-1, 0, 0),vec3(0, 1, 0),
	vec3(0, -1, 0),vec3(0, 0, 1),vec3(0, 0, -1)
};

/***
 *               ___     __        __     __
 *     |\/|  /\   |     |__)  /\  /__` | /  `
 *     |  | /~~\  |     |__) /~~\ .__/ | \__,
 *
 */

inline bool is_identity_auvp(const mat<3, 3>& m, double eps = 1e-15) {
	// [BL] TODO: test extra-diagonal elements ? (can be a big cow without that !!!)
	return std::abs(m[0][0] - 1.) < eps && std::abs(m[1][1] - 1.) < eps && std::abs(m[2][2] - 1.) < eps;
}

mat<3, 3> normalize_columns(const mat<3, 3>& B);
mat<3, 3> invert_columns_norm(const mat<3, 3>& B);

/***
 *               ___     ___            ___  __      __   __  ___
 *     |\/|  /\   |     |__  |  | |    |__  |__)    |__) /  \  |
 *     |  | /~~\  |     |___ \__/ |___ |___ |  \    |  \ \__/  |
 *
 */

mat<3, 3> rotx(double angle);
mat<3, 3> roty(double angle);
mat<3, 3> rotz(double angle);

// non optimized version is  "return rotz(xyz[2]) *roty(xyz[1]) *rotx(xyz[0]);"
mat<3, 3> euler_to_mat(vec3 xyz);

//http://www.staff.city.ac.uk/~sbbh653/publications/euler.pdf    
vec3 mat_to_euler(const mat<3, 3>& r);

/***
 *                __      __   ___  __             ___      ___    __
 *     /\  \_/ | /__`    |__) |__  |__)  |\/| |  |  |   /\   |  | /  \ |\ |
 *    /~~\ / \ | .__/    |    |___ |  \  |  | \__/  |  /~~\  |  | \__/ | \|
 *
 */


struct AxisPermutation {
	AxisPermutation(int id = 0) { mid = id; }
	void aligns_B_wrt_ref(mat<3, 3> ref, mat<3, 3> B);
	void make_col2_equal_to_z(mat<3, 3> B, vec3 z);
	const mat<3, 3>& get_mat() const;
	bool is_identity() { return mid == 0; }
	double operator()(int i, int j) { return get_mat()[i][j]; }
    int stable_axis();
	AxisPermutation inverse();

	int mid;
};

AxisPermutation operator*(AxisPermutation A, AxisPermutation B);
inline vec3 operator*(const AxisPermutation& p, const vec3& v) { return p.get_mat() * v; }

/***
 *     ___  __              ___
 *    |__  |__)  /\   |\/| |__
 *    |    |  \ /~~\  |  | |___
 *
 */

struct Frame {
	mat<3, 3>& r;
	Frame(mat<3, 3>& M) :r(M) {}

	inline mat<3, 3> apply_permutation(AxisPermutation& ap) { return r * ap.get_mat(); }

	void make_z_equal_to(vec3 z);

	static mat<3, 3> average_frame(std::vector<mat<3, 3>>& data);
	static mat<3, 3> representative_frame(std::vector<vec3>& bunch_of_vectors, std::vector<double>& w);
	static mat<3, 3> representative_frame(std::vector<vec3>& bunch_of_vectors);
};


