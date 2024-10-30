#ifndef GEOM2__H__
#define  GEOM2__H__

#include <ultimaille/all.h>
#include <basic.h>

using namespace UM;
#undef max

namespace Geom2 {

	static const mat2x2 Rot90_inv[4] = {
		{vec2(1,0),vec2(0,1)},
		{vec2(0,1),vec2(-1,0)},
		{vec2(-1,0),vec2(0,-1)},
		{vec2(0,-1),vec2(1,0)}
	};
	static const mat2x2 Rot90[4] = {
		{vec2(1,0),vec2(0,1)},
		{vec2(0,-1),vec2(1,0)},
		{vec2(-1,0),vec2(0,-1)},
		{vec2(0,1),vec2(-1,0)}
	};

	static inline mat2x2 rot(double angle) {
		double c = std::cos(angle);
		double s = std::sin(angle);
		return { vec2(c,-s),vec2(s,c) };
	}
};


#endif