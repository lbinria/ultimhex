
#include "frame3D.h"
#include <basic.h>
using namespace UM;


Quaternion axis_rot_to_quaternion(const vec3 axis, const double alpha) {
	Quaternion q; q.v = axis * sin(alpha * .5); q.w = cos(alpha * .5); return q;
}

Quaternion align_with_uv(vec3 u, vec3 v) {
	Quaternion q;
	v.normalize();
	u.normalize();
	vec3 inbetwen_uv(v + u);
	if (inbetwen_uv.norm2() < 1e-10) {
		q.v = vec3(0, 0, 0);
		q.w = -1;
		return q;
	}
	inbetwen_uv.normalize();
	q.w = v * inbetwen_uv; // scalar product with (1,0,0) divided by norm
	q.v =  cross(inbetwen_uv, v); // cross product with (1,0,0) 
	return q;
}

inline mat<3, 3> mat_from_coeffs(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22) {
	mat<3, 3> res;
	res.rows[0] = vec3(a00, a01, a02);
	res.rows[1] = vec3(a10, a11, a12);
	res.rows[2] = vec3(a20, a21, a22);
	return res;
}


Quaternion align_with_v(Quaternion const& q, vec3 v, int& axis) {
	mat<3, 3> rot = q.rotation_matrix();
	v.normalize();
	axis = 0;
	double dot_v = -10e10;
	FOR(l, 6) {
		double dot_tmp = v * (rot * label_to_normal[l]);// transform_vector(label_to_normal[l], rot);
		if (dot_tmp > dot_v) { dot_v = dot_tmp; axis = int(l); }
	}
	vec3 ave = (v + label_to_normal[axis]).normalized();
	return axis_rot_to_quaternion(cross(ave, label_to_normal[axis]), label_to_normal[axis] * ave);
}


vec3 quat_to_vec3(Quaternion const& q) {
	double angle = atan2(q.v.norm(), q.w);
	return 2. * angle * q.v / q.v.norm();
}
Quaternion vec3_to_quat(vec3 const& v) {
	double angle = v.norm() ;
	if (angle < 10e-5) return Quaternion();
	vec3 v_unit = v;
	v_unit.normalize();
	return axis_rot_to_quaternion(v_unit,angle);
}


mat<3, 3> normalize_columns(const mat<3, 3>& B) {
	mat<3, 3> res;
	vec3 n;
	FOR(j, 3) n[j] = B.col(j).norm();// col(B, j).length(); //= std::sqrt(pow(B(0, j), 2) + pow(B(1, j), 2) + pow(B(2, j), 2));
	FOR(j, 3) um_assert(n[j] > 1e-20);
	FOR(i, 3)FOR(j, 3) res[i][j] = B[i][j] / n[j];
	return res;
}

mat<3, 3> invert_columns_norm(const mat<3, 3>& B) {
	mat<3, 3> res;
	vec3 n;
	FOR(j, 3) n[j] = B.col(j).norm2();// pow(B(0, j), 2) + pow(B(1, j), 2) + pow(B(2, j), 2);
	FOR(j, 3) um_assert(n[j] > 1e-20);
	FOR(i, 3)FOR(j, 3) res[i][j] = B[i][j] / n[j];
	return res;
}

/**************************************************************************************************/

mat<3, 3> rotx(double angle) {
	double c = cos(angle);
	double s = sin(angle);
	mat<3, 3> res = mat_from_coeffs(
		1, 0, 0,
		0, c, -s,
		0, s, c
	);
	return res;
}

mat<3, 3> roty(double angle) {
	double c = cos(angle);
	double s = sin(angle);
	mat<3, 3> res = mat_from_coeffs(
		c, 0, s,
		0, 1, 0,
		-s, 0, c
	);
	return res;
}

mat<3, 3> rotz(double angle) {
	double c = cos(angle);
	double s = sin(angle);
	mat<3, 3> res = mat_from_coeffs(
		c, -s, 0,
		s, c, 0,
		0, 0, 1
	);
	return res;
}

// non optimized version is  "return rotz(xyz[2]) *roty(xyz[1]) *rotx(xyz[0]);"
mat<3, 3> euler_to_mat(vec3 xyz) {
	double ca = cos(xyz[0]), sa = sin(xyz[0]);
	double cb = cos(xyz[1]), sb = sin(xyz[1]);
	double cg = cos(xyz[2]), sg = sin(xyz[2]);
	mat<3, 3> res = mat_from_coeffs(
		cb * cg, cg * sa * sb - ca * sg, ca * cg * sb + sa * sg,
		cb * sg, sa * sb * sg + ca * cg, ca * sb * sg - cg * sa,
		-sb, cb * sa, ca * cb
	);
	return res;
}

vec3 mat_to_euler(const mat<3, 3>& r) {//http://www.staff.city.ac.uk/~sbbh653/publications/euler.pdf
	vec3 res;
	if (std::abs(std::abs(r[2][0]) - 1) > 1e-5) {
		res[1] = -asin(r[2][0]);
		res[0] = atan2(r[2][1], r[2][2]);
		res[2] = atan2(r[1][0], r[0][0]);
	}
	else {
		res[2] = 0;
		if (std::abs(r[2][0] + 1) < 1e-5) {
			res[1] = M_PI / 2.;
			res[0] = atan2(r[0][1], r[0][2]);
		}
		else {
			res[1] = -M_PI / 2.;
			res[0] = atan2(-r[0][1], -r[0][2]);
		}
	}
	return res;
}


/***********************************************************************************************************************************/


static mat<3, 3> AxisPermutations[24] = {
mat_from_coeffs(1,0,0,0,1,0,0,0,1),
	mat_from_coeffs(0,-1,0,1,0,0,0,0,1),mat_from_coeffs(-1,0,0,0,-1,0,0,0,1),mat_from_coeffs(0,1,0,-1,0,0,0,0,1),mat_from_coeffs(0,0,1,1,0,0,0,1,0),
	mat_from_coeffs(-1,0,0,0,0,1,0,1,0),mat_from_coeffs(0,0,-1,-1,0,0,0,1,0),mat_from_coeffs(1,0,0,0,0,-1,0,1,0),mat_from_coeffs(0,1,0,0,0,1,1,0,0),
	mat_from_coeffs(0,0,-1,0,1,0,1,0,0),mat_from_coeffs(0,-1,0,0,0,-1,1,0,0),mat_from_coeffs(0,0,1,0,-1,0,1,0,0),mat_from_coeffs(0,0,-1,0,-1,0,-1,0,0),
	mat_from_coeffs(0,1,0,0,0,-1,-1,0,0),mat_from_coeffs(0,0,1,0,1,0,-1,0,0),mat_from_coeffs(0,-1,0,0,0,1,-1,0,0),mat_from_coeffs(0,-1,0,-1,0,0,0,0,-1),
	mat_from_coeffs(1,0,0,0,-1,0,0,0,-1),mat_from_coeffs(0,1,0,1,0,0,0,0,-1),mat_from_coeffs(-1,0,0,0,1,0,0,0,-1),mat_from_coeffs(-1,0,0,0,0,-1,0,-1,0),
	mat_from_coeffs(0,0,1,-1,0,0,0,-1,0),mat_from_coeffs(1,0,0,0,0,1,0,-1,0),mat_from_coeffs(0,0,-1,1,0,0,0,-1,0)
};

static int AxisPermutations_inv[24] = { 0,3,2,1,8,5,15,22,4,14,21,11,12,23,9,6,16,17,18,19,20,10,7,13 };

static int AxisPermutation_mult[24][24] = {
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23},
{1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,21,22,23,20},
{2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,18,19,16,17,22,23,20,21},
{3,0,1,2,7,4,5,6,11,8,9,10,15,12,13,14,19,16,17,18,23,20,21,22},
{4,11,21,14,8,3,13,18,0,7,17,22,20,19,5,2,12,23,9,6,16,15,1,10},
{5,8,22,15,9,0,14,19,1,4,18,23,21,16,6,3,13,20,10,7,17,12,2,11},
{6,9,23,12,10,1,15,16,2,5,19,20,22,17,7,0,14,21,11,4,18,13,3,8},
{7,10,20,13,11,2,12,17,3,6,16,21,23,18,4,1,15,22,8,5,19,14,0,9},
{8,22,15,5,0,14,19,9,4,18,23,1,16,6,3,21,20,10,7,13,12,2,11,17},
{9,23,12,6,1,15,16,10,5,19,20,2,17,7,0,22,21,11,4,14,13,3,8,18},
{10,20,13,7,2,12,17,11,6,16,21,3,18,4,1,23,22,8,5,15,14,0,9,19},
{11,21,14,4,3,13,18,8,7,17,22,0,19,5,2,20,23,9,6,12,15,1,10,16},
{12,6,9,23,16,10,1,15,20,2,5,19,0,22,17,7,4,14,21,11,8,18,13,3},
{13,7,10,20,17,11,2,12,21,3,6,16,1,23,18,4,5,15,22,8,9,19,14,0},
{14,4,11,21,18,8,3,13,22,0,7,17,2,20,19,5,6,12,23,9,10,16,15,1},
{15,5,8,22,19,9,0,14,23,1,4,18,3,21,16,6,7,13,20,10,11,17,12,2},
{16,19,18,17,20,23,22,21,12,15,14,13,8,11,10,9,0,3,2,1,4,7,6,5},
{17,16,19,18,21,20,23,22,13,12,15,14,9,8,11,10,1,0,3,2,5,4,7,6},
{18,17,16,19,22,21,20,23,14,13,12,15,10,9,8,11,2,1,0,3,6,5,4,7},
{19,18,17,16,23,22,21,20,15,14,13,12,11,10,9,8,3,2,1,0,7,6,5,4},
{20,13,7,10,12,17,11,2,16,21,3,6,4,1,23,18,8,5,15,22,0,9,19,14},
{21,14,4,11,13,18,8,3,17,22,0,7,5,2,20,19,9,6,12,23,1,10,16,15},
{22,15,5,8,14,19,9,0,18,23,1,4,6,3,21,16,10,7,13,20,2,11,17,12},
{23,12,6,9,15,16,10,1,19,20,2,5,7,0,22,17,11,4,14,21,3,8,18,13}
};
static int AxisPermutation_stable_direction[24] = {
2,2,2,2,-1,-1,-1,0,-1,1,-1,-1,-1,-1,1,-1,-1,0,-1,1,-1,-1,0,-1 };

//void generate_code() {
//	std::ofstream myfile;
//	myfile.open("C:/NICO/code.txt");
//	myfile << "static int AxisPermutation_mult[24][24] = {\n";
//	FOR(i, 24) {
//		myfile << "{";
//		FOR(j, 24) {
//			mat<3, 3> rot = AxisPermutation(i).get_mat() * AxisPermutation(j).get_mat();
//			FOR(res, 24) if ((rot - AxisPermutation(res).get_mat()).norm() == 0) myfile << res;
//			if (j < 23)		myfile << ",";
//		}
//		if (i < 23)		myfile << "},\n";
//		else myfile << "}\n";
//	}
//	myfile << "};\n";
//
//	myfile << "static int AxisPermutation_stable_direction[24] = {\n";
//	FOR(i, 24) {
//		mat<3, 3> rot = AxisPermutation(i).get_mat();
//		int res = -1;
//		FOR(d, 3) if ((rot.col(d) - mat<3, 3>::identity().col(d)).norm2() == 0) res = d;
//		myfile << res;
//		if (i < 23)		myfile << ",";
//	}
//	myfile << "};\n";
//	myfile.close();
//}


AxisPermutation operator*(AxisPermutation A, AxisPermutation B) {
	return AxisPermutation(AxisPermutation_mult[A.mid][B.mid]);
}
int AxisPermutation::stable_axis() {return AxisPermutation_stable_direction[mid];}

void AxisPermutation::aligns_B_wrt_ref(mat<3, 3> ref, mat<3, 3> B) {
	ref = normalize_columns(ref);
	B = normalize_columns(B);
	int best_i = 0;
	double best_score = 1e20;
	FOR(i, 24) {
		mat<3, 3> m = B * AxisPermutations[i] - ref;
		double score = m.rows[0] * m.rows[0] + m.rows[1] * m.rows[1] + m.rows[2] * m.rows[2];
		if (score < best_score) {
			best_i = i;
			best_score = score;
		}
	}
	mid = best_i;
}

void AxisPermutation::make_col2_equal_to_z(mat<3, 3> B, vec3 z) {
	B = normalize_columns(B);
	int best_i = 0;
	double best_score = 1e20;
	FOR(i, 24) {
		double score = ((B * AxisPermutations[i]).col(2) - z).norm2();
		if (score < best_score) {
			best_i = i;
			best_score = score;
		}
	}
	mid = best_i;
}

const mat<3, 3>& AxisPermutation::get_mat() const {
	return AxisPermutations[mid];
}

AxisPermutation AxisPermutation::inverse() {
	return AxisPermutation(AxisPermutations_inv[mid]);
}

/***********************************************************************************************************************************/

void Frame::make_z_equal_to(vec3 z) {
	z.normalize();
	vec3 x(1, 0, 0);
	vec3 y = cross(z, x);
	if (y.norm2() < .1) {
		x = vec3(0, 1, 0);
		y = cross(z, x);
	}
	y.normalize();
	x = cross(y, z);
	FOR(d, 3) r[d][0] = x[d];
	FOR(d, 3) r[d][1] = y[d];
	FOR(d, 3) r[d][2] = z[d];
}

mat<3, 3> Frame::average_frame(std::vector<mat<3, 3>>& data) {
	SH4 sum;
	FOR(i, data.size()) {
		SH4 nv;
		nv[4] = std::sqrt(7. / 12.);
		nv[8] = std::sqrt(5. / 12.);
		nv.euler_rot(mat_to_euler(data[i]));
		sum = sum + nv;
	}
	sum = sum * (1. / sum.norm());
	return sum.project_mat();
}

mat<3, 3> Frame::representative_frame(std::vector<vec3>& bunch_of_vectors, std::vector<double>& w) {
	SH4 sum;
	FOR(i, bunch_of_vectors.size()) {
		vec3 n = bunch_of_vectors[i];
		n.normalize();
		mat<3, 3> r;
		Frame(r).make_z_equal_to(n);
		SH4 nv;
		nv[4] = std::sqrt(7. / 12.);
		nv.euler_rot(mat_to_euler(r));
		nv = nv * w[i];
		sum = sum + nv;
	}
	sum = sum * (1. / sum.norm());
	return sum.project_mat();
}

mat<3, 3> Frame::representative_frame(std::vector<vec3>& bunch_of_vectors) {
	std::vector<double> w(bunch_of_vectors.size(), 1);
	return  representative_frame(bunch_of_vectors, w);
}



