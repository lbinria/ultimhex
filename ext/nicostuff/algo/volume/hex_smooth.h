#ifndef HEX_SMOOTH__H__
#define	HEX_SMOOTH__H__

#include <ultimaille/all.h>
#include <framework/trace.h>

#include <iostream>


#include <algo/volume/SDF.h>







struct HexSmoother {
	struct TangentConstraint {
		int v;
		vec3 P;
		vec3 n;
	};

	HexSmoother(Hexahedra& m):m(m),constraint_weight(0), lock(m){}
	void set_lock(PointAttribute<bool>& locked_vertex,bool value = true) { FOR(v, m.nverts()) lock[v] = (value==locked_vertex[v]); }
	void add_constraint(int v, vec3 P, vec3 n) { tan_constraints.push_back({v,P,n}); }

	void smooth_LS(double constraint_weight = 100.);
	void smooth_elliptic(double constraint_weight = 100.);

	void show_scaled_jacobien(std::string name, double max_value = 1.1);

	Hexahedra& m;
	PointAttribute<bool> lock;
	std::vector<TangentConstraint> tan_constraints;
	double constraint_weight ;
};

#endif
