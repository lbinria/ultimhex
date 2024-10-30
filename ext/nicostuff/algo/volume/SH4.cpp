
#include <cmath>

#include "SH4.h"
#include "frame3D.h"
#include <math.h>


//    ____            _                     _                  _   _   _                                              _          _       _  _
//   / ___|   _ __   | |__     ___   _ __  (_)   ___    __ _  | | | | | |   __ _   _ __   _ __ ___     ___    _ __   (_)   ___  | |     | || |
//   \___ \  | '_ \  | '_ \   / _ \ | '__| | |  / __|  / _` | | | | |_| |  / _` | | '__| | '_ ` _ \   / _ \  | '_ \  | |  / __| | |     | || |_
//    ___) | | |_) | | | | | |  __/ | |    | | | (__  | (_| | | | |  _  | | (_| | | |    | | | | | | | (_) | | | | | | | | (__  | |___  |__   _|
//   |____/  | .__/  |_| |_|  \___| |_|    |_|  \___|  \__,_| |_| |_| |_|  \__,_| |_|    |_| |_| |_|  \___/  |_| |_| |_|  \___| |_____|    |_|
//           |_|

double SH4::basis(int id, const vec3& v) {
	double x = v.x, x2 = x * x;
	double y = v.y, y2 = y * y;
	double z = v.z, z2 = z * z;
	if (id == 0) return (3. / 4.) * std::sqrt(35. / M_PI) * x * y * (x2 - y2);
	if (id == 1) return (3. / 4.) * std::sqrt(35. / (2. * M_PI)) * z * y * (3 * x2 - y2);
	if (id == 2) return (3. / 4.) * std::sqrt(5. / M_PI) * x * y * (7. * z2 - 1);
	if (id == 3) return (3. / 4.) * std::sqrt(5. / (2. * M_PI)) * z * y * (7. * z2 - 3);
	if (id == 4) return (3. / 16.) * std::sqrt(1. / M_PI) * (35. * z2 * z2 - 30. * z2 + 3);
	if (id == 5) return (3. / 4.) * std::sqrt(5. / (2. * M_PI)) * z * x * (7. * z * z - 3);
	if (id == 6) return (3. / 8.) * std::sqrt(5. / M_PI) * (x2 - y2) * (7. * z * z - 1);
	if (id == 7) return (3. / 4.) * std::sqrt(35. / (2. * M_PI)) * z * x * (x2 - 3. * y2);
	// id==8
	return (3. / 16.) * std::sqrt(35. / M_PI) * (x2 * (x2 - 3. * y2) - y2 * (3. * x2 - y2));
}

void SH4::Rz(double a) {
	SH4 copy(*this);
	double s, c;
	FOR(i, 4) {
		s = sin(a * (4 - i));
		c = cos(a * (4 - i));
		coeff[i] = copy[8 - i] * s + copy[i] * c;
		coeff[8 - i] = copy[8 - i] * c - copy[i] * s;
	}
}

void SH4::Ry(double alpha) {
	SH4 c(*this);
	double sa = sin(alpha), ca = cos(alpha);
	double s2a = sin(2. * alpha), c2a = cos(2. * alpha);
	double s3a = sin(3. * alpha), c3a = cos(3. * alpha);
	double s4a = sin(4. * alpha), c4a = cos(4. * alpha);
	coeff[0] = (c3a + ca * 7.) * .125 * c[0] + (3. * s3a + 7. * sa) * sqrt(.0078125) * c[1] + -(c3a - ca) * sqrt(.109375) * c[2] + -(s3a - 3. * sa) * sqrt(.0546875) * c[3];
	coeff[1] = -(3. * s3a + 7. * sa) * sqrt(.0078125) * c[0] + (9. * c3a + 7. * ca) * .0625 * c[1] + (3. * s3a - sa) * sqrt(.0546875) * c[2] + -(c3a - ca) * sqrt(.24609375) * c[3];
	coeff[2] = -(c3a - ca) * sqrt(.109375) * c[0] + -(3. * s3a - sa) * sqrt(.0546875) * c[1] + (7. * c3a + ca) * .125 * c[2] + (7. * s3a + 3. * sa) * sqrt(.0078125) * c[3];
	coeff[3] = (s3a - 3. * sa) * sqrt(.0546875) * c[0] + -(c3a - ca) * sqrt(.24609375) * c[1] + -(7. * s3a + 3. * sa) * sqrt(.0078125) * c[2] + (7. * c3a + 9. * ca) * .0625 * c[3];
	coeff[4] = (35. * c4a + 20. * c2a + 9.) * .015625 * c[4] + -(7. * s4a + 2. * s2a) * sqrt(.009765625) * c[5] + -(7. * c4a - 4. * c2a - 3.) * sqrt(.0048828125) * c[6] + (s4a - 2. * s2a) * sqrt(.068359375) * c[7] + (c4a - 4. * c2a + 3.) * sqrt(.008544921875) * c[8];
	coeff[5] = (7. * s4a + 2. * s2a) * sqrt(.009765625) * c[4] + (7. * c4a + c2a) * .125 * c[5] + -(7. * s4a - 2. * s2a) * sqrt(.0078125) * c[6] + -(c4a - c2a) * sqrt(.109375) * c[7] + (s4a - 2. * s2a) * sqrt(.013671875) * c[8];
	coeff[6] = -(7. * c4a - 4. * c2a - 3.) * sqrt(.0048828125) * c[4] + (7. * s4a - 2. * s2a) * sqrt(.0078125) * c[5] + (7. * c4a + 4. * c2a + 5.) * .0625 * c[6] + -(s4a + 2. * s2a) * sqrt(.0546875) * c[7] + -(c4a + 4. * c2a - 5.) * sqrt(.0068359375) * c[8];
	coeff[7] = -(s4a - 2. * s2a) * sqrt(.068359375) * c[4] + -(c4a - c2a) * sqrt(.109375) * c[5] + (s4a + 2. * s2a) * sqrt(.0546875) * c[6] + (c4a + 7. * c2a) * .125 * c[7] + -(s4a + 14. * s2a) * sqrt(.001953125) * c[8];
	coeff[8] = (c4a - 4. * c2a + 3.) * sqrt(.008544921875) * c[4] + -(s4a - 2. * s2a) * sqrt(.013671875) * c[5] + -(c4a + 4. * c2a - 5.) * sqrt(.0068359375) * c[6] + (s4a + 14. * s2a) * sqrt(.001953125) * c[7] + (c4a + 28. * c2a + 35.) * .015625 * c[8];
}

void SH4::Rx(double alpha) {
	SH4 c(*this);
	double sa = sin(alpha), ca = cos(alpha);
	double s2a = sin(2. * alpha), c2a = cos(2. * alpha);
	double s3a = sin(3. * alpha), c3a = cos(3. * alpha);
	double s4a = sin(4. * alpha), c4a = cos(4. * alpha);
	coeff[0] = (c3a + ca * 7.) * .125 * c[0] + (c3a - ca) * sqrt(.109375) * c[2] + -(s3a - 3. * sa) * sqrt(.0546875) * c[5] + -(3. * s3a + 7. * sa) * sqrt(.0078125) * c[7];
	coeff[1] = (c4a + 7. * c2a) * .125 * c[1] + (c4a - c2a) * sqrt(.109375) * c[3] + -(s4a - 2. * s2a) * sqrt(.068359375) * c[4] + -(s4a + 2. * s2a) * sqrt(.0546875) * c[6] + -(s4a + 14. * s2a) * sqrt(.001953125) * c[8];
	coeff[2] = (c3a - ca) * sqrt(.109375) * c[0] + (7. * c3a + ca) * .125 * c[2] + -(7. * s3a + 3. * sa) * sqrt(.0078125) * c[5] + -(3. * s3a - sa) * sqrt(.0546875) * c[7];
	coeff[3] = (c4a - c2a) * sqrt(.109375) * c[1] + (7. * c4a + c2a) * .125 * c[3] + -(7. * s4a + 2. * s2a) * sqrt(.009765625) * c[4] + -(7. * s4a - 2. * s2a) * sqrt(.0078125) * c[6] + -(s4a - 2. * s2a) * sqrt(.013671875) * c[8];
	coeff[4] = (s4a - 2. * s2a) * sqrt(.068359375) * c[1] + (7. * s4a + 2. * s2a) * sqrt(.009765625) * c[3] + (35. * c4a + 20. * c2a + 9.) * .015625 * c[4] + (7. * c4a - 4. * c2a - 3.) * sqrt(.0048828125) * c[6] + (c4a - 4. * c2a + 3.) * sqrt(.008544921875) * c[8];
	coeff[5] = (s3a - 3. * sa) * sqrt(.0546875) * c[0] + (7. * s3a + 3 * sa) * sqrt(.0078125) * c[2] + (7. * c3a + 9. * ca) * .0625 * c[5] + (c3a - ca) * sqrt(.24609375) * c[7];
	coeff[6] = (s4a + 2. * s2a) * sqrt(.0546875) * c[1] + (7. * s4a - 2. * s2a) * sqrt(.0078125) * c[3] + (7. * c4a - 4. * c2a - 3.) * sqrt(.0048828125) * c[4] + (7. * c4a + 4. * c2a + 5.) * .0625 * c[6] + (c4a + 4. * c2a - 5.) * sqrt(.0068359375) * c[8];
	coeff[7] = (3. * s3a + 7. * sa) * sqrt(.0078125) * c[0] + (3 * s3a - sa) * sqrt(.0546875) * c[2] + (c3a - ca) * sqrt(.24609375) * c[5] + (9. * c3a + 7. * ca) * .0625 * c[7];
	coeff[8] = (s4a + 14. * s2a) * sqrt(.001953125) * c[1] + (s4a - 2. * s2a) * sqrt(.013671875) * c[3] + (c4a - 4. * c2a + 3.) * sqrt(.008544921875) * c[4] + (c4a + 4. * c2a - 5.) * sqrt(.0068359375) * c[6] + (c4a + 28. * c2a + 35.) * .015625 * c[8];
}

SH4 SH4::Ex() const {
	return SH4(-sqrt(2.) * coeff[7], -sqrt(2.) * coeff[8] - sqrt(3.5) * coeff[6], -sqrt(3.5) * coeff[7] - sqrt(4.5) * coeff[5], -sqrt(4.5) * coeff[6] - sqrt(10.) * coeff[4], sqrt(10.) * coeff[3], sqrt(4.5) * coeff[2], sqrt(3.5) * coeff[1] + sqrt(4.5) * coeff[3], sqrt(2.) * coeff[0] + sqrt(3.5) * coeff[2], sqrt(2.) * coeff[1]);
}

SH4 SH4::Ey() const {
	return SH4(sqrt(2.) * coeff[1], -sqrt(2.) * coeff[0] + sqrt(3.5) * coeff[2], -sqrt(3.5) * coeff[1] + sqrt(4.5) * coeff[3], -sqrt(4.5) * coeff[2], -sqrt(10.) * coeff[5], -sqrt(4.5) * coeff[6] + sqrt(10.) * coeff[4], -sqrt(3.5) * coeff[7] + sqrt(4.5) * coeff[5], -sqrt(2.) * coeff[8] + sqrt(3.5) * coeff[6], sqrt(2.) * coeff[7]);
}

SH4 SH4::Ez() const {
	return SH4(4 * coeff[8], 3 * coeff[7], 2 * coeff[6], coeff[5], 0, -coeff[3], -2 * coeff[2], -3 * coeff[1], -4 * coeff[0]);
}


mat<3, 3> SH4::project_mat(double grad_threshold, double dot_threshold, vec3* euler_prev) {
	SH4 init_harmonics[5] = { rest_frame(),rest_frame(),rest_frame(),rest_frame(),rest_frame() };
	vec3 init_rot[5] = { vec3(0, 0, 0),vec3(M_PI / 4., 0, 0),vec3(0, M_PI / 4., 0),vec3(0, 0, M_PI / 4.),vec3(M_PI / 4., 0, M_PI / 4.) };
	FOR(i, 5) init_harmonics[i].euler_rot(init_rot[i]);

	mat<3, 3> W;
	SH4 v;
	double dot = -1.;

	SH4 query = *this;
	query = query / query.norm();


	if (euler_prev) {
		SH4 prev_seed = init_harmonics[0];

		prev_seed.euler_rot(*euler_prev);
		double tdot = prev_seed * query;
		if (tdot > dot) {
			dot = tdot;

			W = euler_to_mat(*euler_prev);
			v = prev_seed;
		}
	}


	FOR(i, 5) {
		double tdot = init_harmonics[i] * query;
		if (tdot > dot) {
			dot = tdot;

			W = euler_to_mat(init_rot[i]);
			v = init_harmonics[i];
		}
	}

	int cnt = 0;
	double olddot = dot;
	while (cnt < 10000) {
		vec3 grad = vec3(query * v.Ex(), query * v.Ey(), query * v.Ez());
		if (grad.norm() < grad_threshold) break;
		grad /= 8.; // constante au pif trouvée de manière experimentale ; en dessous de ça pb de convergence, au-dessus plus lent
		v.Rx(grad[0]); v.Ry(grad[1]); v.Rz(grad[2]);

		W = rotz(grad[2]) * roty(grad[1]) * rotx(grad[0]) * W;
		cnt++;
		dot = v * query;

		if (dot - olddot < dot_threshold) break;
		olddot = dot;
	}
	if (cnt == 10000) std::cerr << "[error] SH projection infinite loop protection" << std::endl;

	return W;
}



