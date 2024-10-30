#ifndef __SURFACE_OPERATORS__H__
#define __SURFACE_OPERATORS__H__
#include <ultimaille/all.h>
#include "toolbox.h"

#include <iostream>



//                                                       _                
//                                                      | |               
//   _ __    ___     __ _   ___   ___   _ __ ___    ___ | |_  _ __  _   _ 
//  | '_ \  / _ \   / _` | / _ \ / _ \ | '_ ` _ \  / _ \| __|| '__|| | | |
//  | | | || (_) | | (_| ||  __/| (_) || | | | | ||  __/| |_ | |   | |_| |
//  |_| |_| \___/   \__, | \___| \___/ |_| |_| |_| \___| \__||_|    \__, |
//                   __/ |                                           __/ |
//                  |___/                                           |___/                       

struct TrianglesOperator {
	TrianglesOperator(Triangles& m) :m(m) { }
	Triangles& m;
	std::map<int, int> facet_new2old;
	template<class T> void transfert_facet_attribute(FacetAttribute<T>& attr) {
		for (auto p : facet_new2old) attr[p.first] = attr[p.second];
	} 
	std::map<int, int> corner_new2old;
	template<class T> void transfert_corner_attribute(CornerAttribute<T>& attr) {
		for (auto p : corner_new2old) attr[p.first] = attr[p.second];
	}
};

struct TrianglesEdgeFlip : TrianglesOperator {
	TrianglesEdgeFlip(Triangles& m, Surface::Halfedge h) :TrianglesOperator(m), h(h), opp(h.opposite()) { }
	bool valid();
	void apply();
	Surface::Halfedge h;
	Surface::Halfedge opp;
};

struct TrianglesEdgeCollapse : TrianglesOperator {
	TrianglesEdgeCollapse(Triangles& m, Surface::Halfedge h) :TrianglesOperator(m), h(h), opp(h.opposite()) { }
	bool valid();
	void apply();
	Surface::Halfedge h;
	Surface::Halfedge opp;
};


//                                                                  _                
//                                                                 | |               
//   _ __ ___    __ _  _ __     __ _   ___   ___   _ __ ___    ___ | |_  _ __  _   _ 
//  | '_ ` _ \  / _` || '_ \   / _` | / _ \ / _ \ | '_ ` _ \  / _ \| __|| '__|| | | |
//  | | | | | || (_| || |_) | | (_| ||  __/| (_) || | | | | ||  __/| |_ | |   | |_| |
//  |_| |_| |_| \__,_|| .__/   \__, | \___| \___/ |_| |_| |_| \___| \__||_|    \__, |
//                    | |       __/ |                                           __/ |
//                    |_|      |___/                                           |___/


struct MapTrianglesEdgeFlip : TrianglesEdgeFlip {
	MapTrianglesEdgeFlip(Triangles& m, Surface::Halfedge h, CornerAttribute<vec2>& U) :TrianglesEdgeFlip(m, h), U(U) { }
	bool valid();
	bool good_for_delaunay() ;
	void apply();
	CornerAttribute<vec2>& U;
};


struct MapTrianglesEdgeCollapse : TrianglesEdgeCollapse {
	MapTrianglesEdgeCollapse(Triangles& m, Surface::Halfedge h, CornerAttribute<vec2>& U) :TrianglesEdgeCollapse(m, h), U(U) { }
	bool valid();
	void apply();
	CornerAttribute<vec2>& U;
};

#endif