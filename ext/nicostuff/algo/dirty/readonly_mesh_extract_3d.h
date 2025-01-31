#include "basic.h"
#include <ultimaille/all.h>
using namespace UM;

struct ReadOnlyMeshExtract3d {

	struct AxisDirection {
		AxisDirection(Volume::Facet f) : id(f % 6)	{}
		AxisDirection(int id) : id(id)					{}
		vec3 vector()									{ return dir_vec[id]; }
		AxisDirection inverse()							{ return AxisDirection(dir_inv[id]); }
		operator int() const { return id; }

		static vec3 dir_vec[6];
		static int dir_inv[6];
		static int cross[6][6];
		int id;
	};

	struct AxisPermutation {
		AxisPermutation(int p_id)							{ id = p_id; }
		mat3x3 matrix()										{ return rot_matrix[id]; }
		AxisPermutation inverse() const						{ return AxisPermutation(transition_rot_inverse[id]); }
		AxisPermutation operator*(AxisPermutation other)	{ return AxisPermutation(transition_rot_times[id][other.id]); }
		vec3 operator*(vec3 v)								{ return matrix() * v; }
		AxisDirection operator*(AxisDirection dir)			{ return AxisDirection(transition_rot_times_direction[id][dir.id]); }
		operator int() const { return id; }
		static mat3x3 rot_matrix[24];
		static int transition_rot_inverse[24];
		static int transition_rot_times[24][24];
		static int transition_rot_times_direction[24][6];
		int id;
	};

	struct TransitionFunction {// f(x) = rot(x) + t
		TransitionFunction(int m_id = 0, vec3 t = vec3(0, 0,0)) : m(m_id), t(t) {}

		// how to transform h.opposite().facet() to glue it
		// if the rotation is undefined, we assume it to be identity 
		TransitionFunction(Volume::Facet f, CellCornerAttribute<vec3>& uvw) : m(0),t(vec3(0,0,0)) {
			Volume::Halfedge h = f.halfedge(0);
			Volume::Halfedge opp = h.opposite_c();
			um_assert(opp.active());
			if (opp.facet() > f) {
				TransitionFunction tf_opp(opp.facet(), uvw);
				*this = tf_opp.inverted();
				return;
			}
			vec3 v0 = uvw[h.to_corner()] - uvw[h.from_corner()];
			vec3 v1 = uvw[h.prev().from_corner()] - uvw[h.from_corner()];
			vec3 oppv0 = uvw[opp.from_corner()] - uvw[opp.next().from_corner()];
			vec3 oppv1 = uvw[opp.prev().from_corner()] - uvw[opp.next().from_corner()];

			//plop(v0 -oppv0);
			//plop(v1 - oppv1);
			if (v0.norm2() > 1e-20 && v1.norm2() > 1e-20)
				FOR(l, 23) if (
					(m * oppv0 - v0).norm2()
					+ (m * oppv1 - v1).norm2()
					 >
					(AxisPermutation::rot_matrix[l + 1] * oppv0 - v0).norm2()
					+(AxisPermutation::rot_matrix[l + 1] * oppv1 - v1).norm2()
					)
					m.id = l + 1;
			t = uvw[h.from_corner()] - m * uvw[opp.to_corner()];
		}

		TransitionFunction inverted()						{ return TransitionFunction(m.inverse(), -(m.inverse().matrix() * t)) ; }
		TransitionFunction apply(TransitionFunction gp) { return TransitionFunction(m * gp.m, m.matrix() * gp.t + t); }
		vec3 apply(vec3 x)									{ return m.matrix() * x + t; }
		void force_integer_translation()					{ FOR(d, 3) t[d] = std::floor(.5 + t[d]); }
		AxisPermutation m;
		vec3 t;
	};

	struct Link {
		Link(int q=-1, AxisPermutation perm = AxisPermutation(-1)) :m(perm), q_id(q) {}
		Link followed_by(Link other) {
			if (q_id == -1 || other.q_id == -1) return Link(-1, -1);
			return Link(other.q_id, other.m * m);
		}
		bool active() { if (q_id != -1) { um_assert(m != -1); return true; } else  return false; }
		AxisPermutation m;
		int q_id; 
	};

	Volume::Halfedge link_opposite(Volume::Halfedge h);
	Volume::Facet link_opposite_facet(Volume::Facet f);

	struct Embedding {
		int f_id;
		vec3 pos;
	};


	ReadOnlyMeshExtract3d(Tetrahedra& tet, CellCornerAttribute<vec3>& uvw);



	// acces to the input mesh
	Tetrahedron uvw_tet(Volume::Cell t);
	//bool is_singu(Surface::Vertex v);


	// manipulate the quad mesh during simplification
	void create_hex(Volume::Cell t, int i, int j, int k);
	Volume::Cell search_hex(Volume::Cell from_tet, vec3 ijk);
	

	bool can_cancel(Volume::Facet f);
	void cancel(Volume::Facet f);



	/*
	* MAIN STEPS
	*/
	void apply(Hexahedra& hex_out, CellFacetAttribute<int>& emb_out);

	//void sanity_check();
	//CRSMatrix reduce_GP_var(PointAttribute<bool>& lock);
	void reduce_GP_var(ConstrainedLeastSquares& cls, PointAttribute<bool>& lock);
	void pre_process();
	void generate_hexes();
	void connect_hexes();
	void propagate_boundary();
	void untangle();
	//void place_vertices();
	void convert_to_hexes(Hexahedra& hexout, CellFacetAttribute<int>& emb_out);


	void smooth_hex_network();
	void show_hex_network(bool show_hex= true, bool show_ribbons = true);
	//void smooth_quad_network(bool translation=true, bool rotation = true, bool scale = true);
	/*
	* LOW LEVEL ACCESS
	*/


	Tetrahedra& tet;
	CellCornerAttribute<vec3>& uvw;
	CellAttribute<int> off_hex;

	Hexahedra hex;
	// facet attributes
	CellAttribute<vec3> hex_center_U;
	CellAttribute<int> tet_cell;				// the triangles facet it lives on 
	CellAttribute<bool> reverted;
	CellAttribute<bool> active;

	// halfedge attributes
	CellFacetAttribute<Link> link;
	CellFacetAttribute<Embedding> emb;
	double ave_edge_size;
};

