#ifndef GP_BASIC_3D__H__
#define GP_BASIC_3D__H__

#include <ultimaille/all.h>
#include <volume/frame3D.h>
#include <toolbox_tetrahedra.h>	
#include <drop_attribute.h>




#define TRACE_OFF(fct){Trace::Section trace_sec(#fct,Trace::DROP_OFF);fct;}
#define TRACE_ON(fct){Trace::Section trace_sec(#fct,Trace::DROP_UNCHANGED);fct;}







	SparseVector express_with_free_variables(ConstrainedLeastSquares& cls,const LinExpr& le);

	// usage:  auto [already_satisfied,impossible] = constraint_status(cls,le);
	std::tuple<bool,bool> constraint_status(ConstrainedLeastSquares& cls,const LinExpr& le);


inline Volume::Halfedge opposite_on_boundary(Volume::Halfedge h) {
	Volume::Halfedge res = h.opposite_f();
	while (true) {
		auto next = res.opposite_c();
		if (!next.active()) return res;
		res = next.opposite_f();
	}
}

struct TetBoundary {
	inline  TetBoundary(Tetrahedra& tet) :tet(tet), tri_facet_(tet), tri(), tet_facet_(tri) {
		tri.points = tet.points;
		for (auto f_tet : tet.iter_facets()) if (!f_tet.opposite().active()) {
			int f_tri = tri.create_facets(1);
			tet_facet_[f_tri] = f_tet;
			tri_facet_[f_tet] = f_tri;
			FOR(lv, 3) tri.vert(f_tri, lv) = f_tet.vertex(lv);
		}
		tri.connect();
	}
	inline Volume::Facet		tet_facet(int f_tri) { return Volume::Facet(tet, tet_facet_[f_tri]); }
	inline Volume::Halfedge	tet_halfedge(int h_tri) { return Volume::Facet(tet, tet_facet_[h_tri / 3]).halfedge(h_tri % 3); }

	inline Triangles::Facet	tri_facet(int f_tet) { return Triangles::Facet(tri, tri_facet_[f_tet]); }
	inline Triangles::Halfedge	tri_halfedge(int h_tet) { return Triangles::Facet(tri, tri_facet_[h_tet / 3]).halfedge(h_tet % 3); }

	Tetrahedra& tet;
	CellFacetAttribute<int> tri_facet_;

	Triangles tri;
	FacetAttribute<int> tet_facet_;
};

struct HexBoundary {
	inline  HexBoundary(Hexahedra& hex,bool duplicate_vertices=false) :hex(hex), quad_facet_(hex), quad(), hex_facet_(quad) {
		if (duplicate_vertices) 
			ToolBox(quad.points).copy_from(hex.points);
		else 
			quad.points = hex.points;
		
		for (auto f_hex : hex.iter_facets()) if (!f_hex.opposite().active()) {
			int f_quad= quad.create_facets(1);
			hex_facet_[f_quad] = f_hex;
			quad_facet_[f_hex] = f_quad;
			FOR(lv, 4) quad.vert(f_quad, lv) = f_hex.vertex(lv);
		}
		quad.connect();
	}
	inline Volume::Facet		hex_facet(int f_quad) { return Volume::Facet(hex, hex_facet_[f_quad]); }
	inline Volume::Halfedge	hex_halfedge(int h_quad) { return Volume::Facet(hex, hex_facet_[h_quad/ 4]).halfedge(h_quad% 4); }

	inline Quads::Facet	quad_facet(int f_hex) { return Quads::Facet(quad, quad_facet_[f_hex]); }
	inline Quads::Halfedge	quad_halfedge(int h_hex) { return Quads::Facet(quad, quad_facet_[h_hex / 4]).halfedge(h_hex % 4); }

	Hexahedra& hex;
	CellFacetAttribute<int> quad_facet_;

	Quads quad;
	FacetAttribute<int> hex_facet_;
};




// jacobian: chaque ligne est le gradiant d'une coord
mat3x3 uvw_to_jacobian(Volume::Cell c, CellCornerAttribute<vec3>& uvw);
mat<4, 3> jacobian_to_uvw(Volume::Cell c, CellAttribute<mat3x3>& J);

void map_to_jacobian(Tetrahedra& m,CellCornerAttribute<vec3>& U_in, CellAttribute<mat3x3>& J_out);
void map_to_jacobian(Tetrahedra& m,CellCornerAttribute<vec3>& U_in, CellAttribute<mat3x3>& J_out, CellAttribute<vec3>& T_out);
void jacobian_to_map(Tetrahedra& m,CellAttribute<mat3x3>& J_in, CellAttribute<vec3>& T_in,CellCornerAttribute<vec3>& U_out);
mat3x3 closest_rotation(mat3x3 J);

AxisPermutation  permute_Jj_to_approx_Ji(mat<3, 3> Ji, mat<3, 3> Jj) ;


struct GPTransitionFunction {// f(x) = rot(x) + t
	GPTransitionFunction();

	// the transition function express the f.opposite().cell() uvw in the f.cell() uvw.
	// the objective is to grow region where GPTransitionFunction makes each local frame compatible with the seed frame
	void align_rot(Volume::Facet f, mat3x3 M_f_cell, mat3x3 M_opp_cell);

	GPTransitionFunction(Volume::Facet f, CellCornerAttribute<vec3>& U);

	GPTransitionFunction inverted();

	GPTransitionFunction apply(GPTransitionFunction gp);

	void show();

	vec3 apply(vec3 x);
	AxisPermutation ap;
	vec3 t;
};

void test_GPTransitionFunction(Tetrahedra& m, CellCornerAttribute<vec3>& U);


inline Tetrahedron uvw_tet(CellCornerAttribute<vec3>& uvw, Volume::Cell c) { return Tetrahedron(uvw[c.corner(0)], uvw[c.corner(1)], uvw[c.corner(2)], uvw[c.corner(3)]); }
inline Triangle3 uvw_tri(CellCornerAttribute<vec3>& uvw, Volume::Facet f) { return Triangle3(uvw[f.corner(0)], uvw[f.corner(1)], uvw[f.corner(2)]); }


inline bool uvw_singu(Volume::Halfedge h,CellCornerAttribute<vec3>& U){
	auto cir = h;
	AxisPermutation perm(0);
	do {
		if (!cir.opposite_f().facet().opposite().active()) return false;
		GPTransitionFunction tf(cir.opposite_f().facet(),U);
		AxisPermutation& ap = tf.ap;
		perm = perm * ap;//[cir.opposite_f().facet()];
		cir = cir.opposite_f().opposite_c();
		if (!cir.active()) return false;
	} while (cir != h);

	return !perm.is_identity();
}


struct DisjointSetWithNeutral : public DisjointSet {
	inline DisjointSetWithNeutral(int n) : DisjointSet(n + 1) {}

	inline void set_neutral(int i) { merge(i, m_ids.size() - 1); }

	inline int get_sets_id(std::vector<int>& id2setid) {
		int last = int(m_ids.size()) - 1;
		id2setid.resize(last); // prepare the correspondance root id => set id
		int nsets = 0;
		for (int i = 0; i < last; i++) {
			if (i != root(i)) continue;
			id2setid[i] = nsets;
			nsets++;
		}
		for (int i = 0; i < last; i++) // fill the rest
			if (m_ids[i] == last) id2setid[i] = -1;
			else id2setid[i] = id2setid[m_ids[i]];
		return nsets;
	}
};



inline void check_is_seamless(Tetrahedra& m , CellCornerAttribute<vec3>& U){
	for(auto f:m.iter_facets()) {
		if (!f.opposite().active()) continue; 
		GPTransitionFunction tf(f,U);
		for(auto h:f.iter_halfedges()){
			double diff = (U[h.from_corner()]-tf.apply(U[h.opposite_c().to_corner()])).norm();
			if (diff>1e-10) {
				Trace::drop_mesh_is_active = true;
				Drop(m,U).apply("Not Seamless");
				drop_triangle(Triangle3(f.vertex(0).pos(),f.vertex(1).pos(),f.vertex(2).pos()),"map not seamless across triangle");
				throw( std::runtime_error("non seamless map"));
			}
	}
}
}


//inline bool is_seamless(Tetrahedra& m , CellCornerAttribute<vec3>& U){
//	for(auto f:m.iter_facets()) {
//		if (!f.opposite().active()) continue; 
//		GPTransitionFunction tf(f,U);
//		for(auto h:f.iter_halfedges()){
//			double diff = (U[h.from_corner()]-tf.apply(U[h.opposite_c().to_corner()])).norm();
//			if (diff>1e-10) return false;
//	}
//	}
//	return true;
//}
inline bool is_GP(Tetrahedra& m , CellCornerAttribute<vec3>& U){
	for(auto f:m.iter_facets()) {
		if (!f.opposite().active()) continue; 
		GPTransitionFunction tf(f,U);
		
		// check seamless
		for(auto h:f.iter_halfedges()){
			double diff = (U[h.from_corner()]-tf.apply(U[h.opposite_c().to_corner()])).norm();
			if (diff>1e-10) return false;
		}
		// check integer translation
		FOR(d,3) if (std::abs(tf.t[d]- std::round(tf.t[d]))>1e-10) return false;
	}
}



	template<class P> void copy_attribute(GenericAttribute<P>& from, GenericAttribute<P>& to) { FOR(i, from.ptr->data.size()) to[i] = from[i]; }




struct SubTetrahedra {
	inline SubTetrahedra(Tetrahedra& m, Tetrahedra& sub_tet, CellAttribute<bool>& selection) : m(m), proxy(sub_tet), ref_cell(proxy), ref_point(proxy.points) {
		ToolBox(proxy).copy_from(m);
		for (auto c : proxy.iter_cells()) ref_cell[c] = c;
		for (auto v : proxy.iter_vertices()) ref_point[v] = v;
		std::vector<bool> to_kill(proxy.ncells());
		for (auto c : proxy.iter_cells()) to_kill[c] = !selection[c];
		proxy.delete_cells(to_kill);
		proxy.delete_isolated_vertices();
		proxy.connect();
	}

	template <class T, class C>
	void to_proxy(C& attr, C& proxy_attr) { um_assert(false); }

	template <class T>
	void to_proxy(CellAttribute<T>& attr, CellAttribute<T>& proxy_attr) {
		for (auto c : proxy.iter_cells()) proxy_attr[c] = attr[ref_cell[c]];
	}
	template <class T>
	void to_proxy(CellCornerAttribute<T>& attr, CellCornerAttribute<T>& proxy_attr) {
		for (auto c : proxy.iter_corners()) proxy_attr[c] = attr[c % 4 + 4 * ref_cell[c / 4]];
	}
	template <class T>
	void to_proxy(CellFacetAttribute<T>& attr, CellFacetAttribute<T>& proxy_attr) {
		for (auto c : proxy.iter_cells()) FOR(lf,4) proxy_attr[4*c+lf] = attr[4*ref_cell[c]+lf];
	}




	template <class T>
	void from_proxy(CellAttribute<T>& attr, CellAttribute<T>& proxy_attr) {
		for (auto c : proxy.iter_cells()) attr[ref_cell[c]] = proxy_attr[c];
	}
	template <class T>
	void from_proxy(CellCornerAttribute<T>& attr, CellCornerAttribute<T>& proxy_attr) {
		for (auto c : proxy.iter_corners()) attr[c % 4 + 4 * ref_cell[c / 4]] = proxy_attr[c];
	}
	template <class T>
	void from_proxy(CellFacetAttribute<T>& attr, CellFacetAttribute<T>& proxy_attr) {
		for (auto c : proxy.iter_cells()) FOR(lf,4) attr[4*ref_cell[c]+lf] = proxy_attr[4*c+lf] ;
	}

	Tetrahedra& m;
	Tetrahedra& proxy;
	CellAttribute<int> ref_cell;
	PointAttribute<int> ref_point;
};


#endif