#include <surface/mesh_edit.h>

int valence(Surface::Vertex&& v) {
	int ret = 0;
	for (Surface::Halfedge cir : v.iter_halfedges()) ret++;
	return ret;
}


//                                                       _                
//                                                      | |               
//   _ __    ___     __ _   ___   ___   _ __ ___    ___ | |_  _ __  _   _ 
//  | '_ \  / _ \   / _` | / _ \ / _ \ | '_ ` _ \  / _ \| __|| '__|| | | |
//  | | | || (_) | | (_| ||  __/| (_) || | | | | ||  __/| |_ | |   | |_| |
//  |_| |_| \___/   \__, | \___| \___/ |_| |_| |_| \___| \__||_|    \__, |
//                   __/ |                                           __/ |
//                  |___/                                           |___/                       


bool TrianglesEdgeFlip::valid() {
	if (!opp.active()) return false;
	if (h.facet().size() != 3) return false;
	if (opp.facet().size() != 3) return false;
	if (valence(h.to()) < 4) return false;
	if (valence(h.from()) < 4) return false;
	for (auto cir_left : h.next().to().iter_halfedges())
		if (cir_left.to() == opp.next().to())
			return false;;
	return true;
}
void TrianglesEdgeFlip::apply() {
	Surface::Halfedge c[2][3] = {
		{ h.next(),h.prev(),opp.prev() },
		{ opp.next(),opp.prev(),h.prev() },
	};
	FOR(f, 2) {
		auto nf = m.Surface::conn->create_facet({ c[f][0].from(), c[f][1].from(), c[f][2].from() });
		facet_new2old[nf] = c[f][0].facet();
		FOR(lh, 3) corner_new2old[nf.halfedge(lh)] = c[f][lh];
		m.Surface::conn->active[c[f][0].facet()] = false;
	}
}

bool TrianglesEdgeCollapse::valid() {
	if (opp.active())	for (auto cir_from : h.from().iter_halfedges()) if (!cir_from.opposite().active()) return false;

	for (auto cir_from : h.from().iter_halfedges())
		for (auto cir_to : h.to().iter_halfedges())
			if (cir_from.to() == cir_to.to())
				if (cir_from.to() != h.next().to())
					if (!opp.active() || cir_from.to() != opp.next().to())
						return false;
	return true;
}

void TrianglesEdgeCollapse::apply() {
	std::vector<int> h_incident;
	for (auto cir_from : h.from().iter_halfedges())
		h_incident.push_back(cir_from);

	for (int it_id : h_incident) {
		Surface::Halfedge it(m, it_id);

		Surface::Halfedge c[3] = { h.next(),it.next(),it.next().next() };
		//if (c[1].from() == h.to() || c[2].from() == h.to()) continue;
		if (c[1].from() == c[0].from()) continue;
		if (c[2].from() == c[0].from()) continue;
		auto nf = m.Surface::conn->create_facet({ c[0].from(), c[1].from(), c[2].from() });
		facet_new2old[nf] = it.facet();
		FOR(lh, 3) corner_new2old[nf.halfedge(lh)] = c[lh];
	}
	for (int it_id : h_incident)
		m.Surface::conn->active[Surface::Halfedge(m, it_id).facet()] = false;
}







//                                                                  _                
//                                                                 | |               
//   _ __ ___    __ _  _ __     __ _   ___   ___   _ __ ___    ___ | |_  _ __  _   _ 
//  | '_ ` _ \  / _` || '_ \   / _` | / _ \ / _ \ | '_ ` _ \  / _ \| __|| '__|| | | |
//  | | | | | || (_| || |_) | | (_| ||  __/| (_) || | | | | ||  __/| |_ | |   | |_| |
//  |_| |_| |_| \__,_|| .__/   \__, | \___| \___/ |_| |_| |_| \___| \__||_|    \__, |
//                    | |       __/ |                                           __/ |
//                    |_|      |___/                                           |___/


bool MapTrianglesEdgeFlip::valid() {
	// check TOPO
	if (!TrianglesEdgeFlip::valid())			return false;
	// check continuity
	if ((U[h] - U[opp.next()]).norm2() > 1e-20) return false;
	if ((U[h.next()] - U[opp]).norm2() > 1e-20) return false;
	// check orientation
	vec2 q[4] = { U[h], U[opp.prev()], U[h.next()], U[h.prev()] };
	FOR(lv, 4)if (mat2x2({ vec2(q[(lv + 1) % 4] - q[lv]),vec2(q[(lv + 2) % 4] - q[lv]) }).det() < 0) return false;
	return true;
}
bool MapTrianglesEdgeFlip::good_for_delaunay() {
	um_assert(opp.active()) ;
	vec2 P = U[opp.prev()];
	vec2 A = U[h] - P;
	vec2 B = U[h.next()] - P;
	vec2 C = U[h.prev()] - P;

	mat3x3 mat = {
		vec3(A[0],A[1],A.norm2()),
		vec3(B[0],B[1],B.norm2()),
		vec3(C[0],C[1],C.norm2())
	};
	return mat.det() > 1e-10;
}
void MapTrianglesEdgeFlip::apply() {
	TrianglesEdgeFlip::apply();
	transfert_corner_attribute(U);
}



bool MapTrianglesEdgeCollapse::valid() {
	if (!TrianglesEdgeCollapse::valid()) return false;
	// check triangles orientation 
	for (auto cir : h.from().iter_halfedges()) {
		if (cir == h) continue;
		if (cir.prev() == h.opposite()) continue;
		if (mat2x2({ U[cir.next()] - U[h.next()],U[cir.prev()] - U[h.next()] }).det() < 0) return false;
	}
	return true;
}

void MapTrianglesEdgeCollapse::apply() {
	TrianglesEdgeCollapse::apply();
	transfert_corner_attribute(U);
}

