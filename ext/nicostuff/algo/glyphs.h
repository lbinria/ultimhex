#ifndef __GLYPHS_H__
#define __GLYPHS_H__

#include "basic.h"


#include <ultimaille/all.h>
using namespace UM;


namespace Glyphs {

	enum ArrowHead {
		NONE, CONE, SPHERE
	};
	struct ArrowStyle {
		ArrowStyle(int p_resolution, double p_diameter = 0, ArrowHead p_arrow_head_style = CONE) {
			resolution = p_resolution; diameter = p_diameter; arrow_head_style = p_arrow_head_style;
		}
		int resolution = 3;
		double diameter = 0;
		ArrowHead arrow_head_style = CONE;
	};

	inline Quaternion align_with_uv(vec3 u, vec3 v) {
		v.normalize();
		u.normalize();
		if (u * v < -.99999) { Quaternion res;  res.v = vec3(1, 0, 0); res.w = 0;  return res; }
		if (std::abs(u * v) > .99999)  return Quaternion();
		vec3 inbetwen_uv(v + u);
		inbetwen_uv.normalize();
		Quaternion res;
		res.w = v * inbetwen_uv; // scalar product with (1,0,0) divided by norm
		res.v = cross(inbetwen_uv, v); // cross product with (1,0,0) 
		return res;
	}






	struct Arrow {
		Arrow(Polygons& m) :m(m), res(3), r(-1) {}
		Arrow& resolution(int reso) { res = reso; return *this; }
		Arrow& radius(double val) { r = val; return *this; }

		Polygons& m;
		int res;
		double r;

		void apply(const vec3& A, const vec3& B) {
			if (r == -1) r = (B - A).norm() / 10.;

			vec3 n = B - A;
			double l = std::max(1e-5, n.norm());
			n = n / l;


			int offv = m.points.create_points(res * 3 + 1);
			FOR(v, res) m.points[offv + v] = r * vec3(cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0);
			FOR(v, res) m.points[offv + res + v] = m.points[offv + v] + vec3(0, 0, std::max(0., l - 2.5 * r));
			FOR(v, res) m.points[offv + 2 * res + v] = 2. * m.points[offv + v] + vec3(0, 0, std::max(0., l - 5 * r));
			m.points[offv + 3 * res] = vec3(0, 0, l);
			int offf = m.create_facets(res, 4);
			FOR(f, res) {
				m.facets[m.offset[offf + f] + 0] = offv + f;
				m.facets[m.offset[offf + f] + 1] = offv + (f + 1) % res;
				m.facets[m.offset[offf + f] + 2] = offv + (f + 1) % res + res;
				m.facets[m.offset[offf + f] + 3] = offv + f + res;
			}
			offf = m.create_facets(res, 3);
			FOR(f, res) {
				m.facets[m.offset[offf + f] + 0] = offv + f + 2 * res;
				m.facets[m.offset[offf + f] + 1] = offv + (f + 1) % res + 2 * res;
				m.facets[m.offset[offf + f] + 2] = offv + 3 * res;
			}

			Quaternion quat = Glyphs::align_with_uv(vec3(0, 0, 1), n);
			auto M = quat.rotation_matrix();
			for (int v = offv; v < offv + 3 * res + 1; v++)
				m.points[v] = M * m.points[v] + A;
		}
		template <class T> void apply(const vec3& A, const vec3& B, FacetAttribute<T>& val, T value) {
			apply(A, B);
			FOR(f, 2 * res) val[m.nfacets() - f - 1] = value;
		}
		template <class T> void apply(const vec3& A, const vec3& B, PointAttribute<T>& val, T value0, T value1) {
			apply(A, B);
			FOR(v, 3 * res + 1) {
				double c = (m.points[m.nverts() - v - 1] - A) * (B - A) / (B - A).norm2();
				val[m.nverts() - v - 1] = (1. - c) * value0 + c * value1;
			}
		}
	};
	struct Tube {
		Tube(Polygons& m) :m(m), res(3), r(-1) {}
		Tube& resolution(int reso) { res = reso; return *this; }
		Tube& radius(double val) { r = val; return *this; }

		Polygons& m;
		int res;
		double r;

		void apply(const vec3& A, const vec3& B) {
			if (r == -1) r = (B - A).norm() / 10.;
			int offv = m.points.create_points(res * 2);
			FOR(v, res) m.points[offv + v] = r * vec3(cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0);
			FOR(v, res) m.points[offv + res + v] = m.points[offv + v] + vec3(0, 0, 1);
			int offf = m.create_facets(res, 4);
			FOR(f, res) {
				m.facets[m.offset[offf + f] + 0] = offv + f;
				m.facets[m.offset[offf + f] + 1] = offv + (f + 1) % res;
				m.facets[m.offset[offf + f] + 2] = offv + (f + 1) % res + res;
				m.facets[m.offset[offf + f] + 3] = offv + f + res;
			}

			vec3 n = B - A;
			double l = std::max(1e-5, n.norm());
			n = n / l;

			for (int v = offv; v < offv + 2 * res; v++)  m.points[v][2] *= l;

			Quaternion quat = Glyphs::align_with_uv(vec3(0, 0, 1), n);
			auto M = quat.rotation_matrix();
			for (int v = offv; v < offv + 2 * res; v++)
				m.points[v] = M * m.points[v] + A;
		}
		template <class T> void apply(const vec3& A, const vec3& B, FacetAttribute<T>& val, T value) {
			apply(A, B);
			FOR(f, res) val[m.nfacets() - f - 1] = value;
		}
		template <class T> void apply(const vec3& A, const vec3& B, PointAttribute<T>& val, T value0, T value1) {
			apply(A, B);
			FOR(v, 2 * res) {
				double c = (m.points[m.nverts() - v - 1] - A) * (B - A) / (B - A).norm2();
				val[m.nverts() - v - 1] = (1. - c) * value0 + c * value1;
			}
		}
	};




	struct Capsule {
		Capsule(Polygons& m) :m(m), res(3), r(-1) {}
		Capsule& resolution(int reso) { res = reso; return *this; }
		Capsule& radius(double val) { r = val; return *this; }

		Polygons& m;
		int res;
		double r;

		void apply(const vec3& A, const vec3& B) {
			if (r == -1) r = (B - A).norm() / 10.;
			res = std::max(res, 4);
			int beta_res = std::max(res / 2, 3);
			int offv = m.points.create_points(res * (beta_res + 1));
			FOR(v, res) FOR(layer, beta_res + 1)
				m.points[offv + layer * res + v] =
				vec3(cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0);
			FOR(layer, beta_res + 1)
				FOR(v, res) {
				double beta = M_PI / 2. + M_PI * double(layer) / double(beta_res);
				m.points[offv + layer * res + v] *= cos(beta);
				m.points[offv + layer * res + v][2] = sin(beta);
			}
			int offf = m.create_facets(res * beta_res, 4);
			FOR(layer, beta_res)FOR(f, res) {
				m.facets[m.offset[offf + layer * res + f] + 0] = offv + layer * res + f;
				m.facets[m.offset[offf + layer * res + f] + 1] = offv + layer * res + (f + 1) % res;
				m.facets[m.offset[offf + layer * res + f] + 2] = offv + layer * res + (f + 1) % res + res;
				m.facets[m.offset[offf + layer * res + f] + 3] = offv + layer * res + f + res;
			}

			vec3 n = B - A;
			double l = std::max(1e-15, n.norm());
			n = n / l;

			for (int v = offv; v < m.nverts(); v++) FOR(d, 3)m.points[v][d] *= r;
			for (int v = offv; v < m.nverts(); v++) if (m.points[v][2] > 0) m.points[v][2] += l;


			Quaternion quat = Glyphs::align_with_uv(vec3(0, 0, 1), n);
			auto M = quat.rotation_matrix();
			for (int v = offv; v < m.nverts(); v++)
				m.points[v] = M * m.points[v] + A;
		}
		template <class T> void apply(const vec3& A, const vec3& B, FacetAttribute<T>& val, T value) {
			apply(A, B);
			int beta_res = std::max(res / 2, 3);
			FOR(f, res * beta_res) val[m.nfacets() - f - 1] = value;
		}
		template <class T> void apply(const vec3& A, const vec3& B, PointAttribute<T>& val, T value0, T value1) {
			apply(A, B);
			int beta_res = std::max(res / 2, 3);
			FOR(v, res * (beta_res + 1)) {
				double c = (m.points[m.nverts() - v - 1] - A) * (B - A) / (B - A).norm2();
				val[m.nverts() - v - 1] = (1. - c) * value0 + c * value1;
			}
		}
	};


	struct Disk {
		Disk(Polygons& m) :m(m), res(3), r(-1) {}
		Disk& resolution(int reso) { res = reso; return *this; }
		Disk& radius(double val) { r = val; return *this; }

		Polygons& m;
		int res;
		double r;

		void apply(const vec3& O, const vec3& normal) {
			int off_f = m.create_facets(res, 3);
			int off_v = m.points.create_points(res + 1);

			m.points[off_v] = vec3(0, 0, 0);
			FOR(i, res) {
				double angle = 2. * M_PI * double(i) / double(res);
				m.points[off_v + i + 1] = r * vec3(cos(angle), sin(angle), 0);
			}

			FOR(i, res) {
				m.vert(off_f + i, 0) = off_v;
				m.vert(off_f + i, 1) = off_v + 1 + i;
				m.vert(off_f + i, 2) = off_v + 1 + (i + 1) % res;
			}
			Quaternion quat = Glyphs::align_with_uv(vec3(0, 0, 1), normal);
			auto M = quat.rotation_matrix();
			FOR(i, res + 1) m.points[off_v + i] = M * m.points[off_v + i] + O;

		}
		template <class T> void apply(const vec3& O, const vec3& normal, FacetAttribute<T>& val, T value) {
			apply(O, normal);
			FOR(f, res) val[m.nfacets() - f - 1] = value;
		}
		template <class T> void apply(const vec3& O, const vec3& normal, PointAttribute<T>& val, T middle_value, T boundary_value) {
			apply(O, normal);
			val[m.nverts() - res - 2] = middle_value;
			FOR(v, res) {
				val[m.nverts() - v - 1] = boundary_value;
			}
		}
	};

	struct Sphere {
		Sphere(Polygons& m) :m(m), res(3), r(-1) {}
		Sphere& resolution(int reso) { res = reso; return *this; }
		Sphere& radius(double val) { r = val; return *this; }

		Polygons& m;
		int res;
		double r;

		void apply(const vec3& O, const vec3& normal) {
			Capsule(m).radius(r).resolution(res).apply(O, O + 1e-10 * normal);
		}
		template <class T> void apply(const vec3& O, const vec3& normal, FacetAttribute<T>& val, T value) {
			Capsule(m).radius(r).resolution(res).apply(O, O + 1e-10 * normal, val, value);
		}
		template <class T> void apply(const vec3& O, const vec3& normal, PointAttribute<T>& val, T value0, T value1) {
			Capsule(m).radius(r).resolution(res).apply(O, O + 1e-10 * normal);
			int beta_res = std::max(res / 2, 3);
			FOR(v, res * (beta_res + 1)) {
				double c = .5 + .5 * (m.points[m.nverts() - v - 1] - O) * (normal) / (r * normal.norm2());
				val[m.nverts() - v - 1] = (1. - c) * value0 + c * value1;
			}

		}
	};




}


#endif