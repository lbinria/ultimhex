#pragma once

#include <functional>
// #include <utility>

#include <ultimaille/all.h>
#include <geogram/mesh/mesh.h>

#include "my_tet_boundary.h"

namespace um_bindings {

	/**
	 * Get GEO cell edge index from UM volume halfedge index
	 */
	static constexpr GEO::index_t cell_e_from_he[24] = {8,7,11,3,1,10,5,9,0,9,4,8,11,6,10,2,3,2,1,0,4,5,6,7};

	/**
	 * Get GEO cell local facet index from UM volume halfedge index
	 */
	static constexpr GEO::index_t cell_lf_from_he[24] = {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5};

	// TODO add GEO cell local facet from UM volume facet !

	/**
	 * Get UM volume local halfedge index in a cell from GEO cell edge and local facet
	 */
	static constexpr int he_from_cell_e_lf(GEO::index_t e, GEO::index_t lf) {
		// Hard-coded matrix (compute direct access to index using index = e + lf * 12)
		constexpr int he_from_elf[72] = {
			-1, -1, -1, 3, -1, -1, -1, 1, 0, -1,
			-1, 2, -1, 4, -1, -1, -1, 6, -1, -1,
			-1, 7, 5, -1, 8, -1, -1, -1, 10, -1,
			-1, -1, 11, 9, -1, -1, -1, -1, 15, -1,
			-1, -1, 13, -1, -1, -1, 14, 12, 19, 18,
			17, 16, -1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 20, 21, 22, 23, -1, -1,
			-1, -1
		};

		return he_from_elf[e + lf * 12];
	}

	static int um_facet_index_from_geo_facet_index(GEO::index_t f, int n_facet_per_cell) {
		// In GEO 8 facets per cell
		int c_idx = f / 8;
		int lf = (f - c_idx * 8);
		return c_idx * n_facet_per_cell + lf;
	}

	static int geo_facet_index_from_um_facet_index(GEO::index_t f, int n_facet_per_cell) {
		// In GEO 8 facets per cell
		int c_idx = f / n_facet_per_cell;
		int lf = (f - c_idx * n_facet_per_cell);
		return c_idx * 8 + lf;
	}

	static GEO::vec3 geo_vec(UM::vec3 v) {
		return GEO::vec3(v.x, v.y, v.z);
	}

	template<typename T>
	void set_um_vertices_from_geo_vertices(GEO::Mesh &m, T &m_out);

	template<typename T>
	void set_geo_vertices_from_um_vertices(T &m, GEO::Mesh &m_out);

	/**
	 * Get ultimaille tri from geogram mesh
	 */
	void um_tri_from_geo_mesh(GEO::Mesh &m, UM::Triangles &m_out);
	
	/**
	 * Get ultimaille tet from geogram mesh
	 */
	void um_tet_from_geo_mesh(GEO::Mesh &m, UM::Tetrahedra &m_out);

	/**
	 * Get ultimaille hex from geogram mesh
	 */
	void um_hex_from_geo_mesh(GEO::Mesh &m, UM::Hexahedra &m_out);

	/**
	 * Get geogram mesh from ultimaille tri
	 */
	void geo_mesh_from_um_tri(UM::Triangles &m, GEO::Mesh &m_out, bool clear = true);

	/**
	 * Get geogram mesh from ultimaille tet
	 */
	void geo_mesh_from_um_tet(UM::Tetrahedra &m, GEO::Mesh &m_out, bool clear = true);

	/**
	 * Get geogram mesh from ultimaille tet + tri surface
	 */
	void geo_mesh_from_tetboundary(TetBoundary &m, GEO::Mesh &m_out);

	/**
	 * Get geogram mesh from ultimaille quad
	 */
	void geo_mesh_from_um_quad(UM::Quads &m, GEO::Mesh &m_out, bool clear);

	/**
	 * Get geogram mesh from ultimaille hex
	 */
	void geo_mesh_from_um_hex(UM::Hexahedra &m, GEO::Mesh &m_out, bool clear = true);

	/**
	 * Get geogram mesh from ultimaille hex + quad surface
	 */
	void geo_mesh_from_hexboundary(HexBoundary &m, GEO::Mesh &m_out);

	/**
	 * Get geogram facet attribute from ultimaille facet attribute
	 */
	template<typename T>
	void geo_attr_from_um_attr(UM::Surface &s, std::string attr_name, UM::FacetAttribute<T> &p, GEO::Mesh &m) {

		using GEO_Attr_T = std::conditional_t<
			std::is_same_v<T, double>, GEO::Numeric::float64,
			std::conditional_t<std::is_same_v<T, int>, GEO::Numeric::uint32,
			std::conditional_t<std::is_same_v<T, bool>, GEO::Numeric::uint8, T>>>;

   		// GEO::Attribute<GEO_Attr_T> geo_attr(
        // 	std::conditional_t<std::is_same_v<T, double>, decltype(m.vertices.attributes()), decltype(m.facets.attributes())>{}(), "flag"
		// );

		GEO::Attribute<GEO_Attr_T> geo_attr(
			m.facets.attributes(), attr_name
		);

		for (auto f : s.iter_facets()) {			
			geo_attr[f] = static_cast<GEO_Attr_T>(p[f]);
		}
	}

	/**
	 * Get geogram facet attribute from ultimaille facet attribute
	 */
	template<typename T>
	void geo_attr_from_um_attr2(UM::Surface &s, std::string attr_name, std::vector<T> &p, GEO::Mesh &m_out) {

		using GEO_Attr_T = std::conditional_t<
			std::is_same_v<T, double>, GEO::Numeric::float64,
			std::conditional_t<std::is_same_v<T, int>, GEO::Numeric::uint32,
			std::conditional_t<std::is_same_v<T, bool>, GEO::Numeric::uint8, T>>>;

   		// GEO::Attribute<GEO_Attr_T> geo_attr(
        // 	std::conditional_t<std::is_same_v<T, double>, decltype(m.vertices.attributes()), decltype(m.facets.attributes())>{}(), "flag"
		// );

		GEO::Attribute<GEO_Attr_T> geo_attr(
			m_out.facets.attributes(), attr_name
		);

		for (auto f : s.iter_facets()) {			
			geo_attr[f] = static_cast<GEO_Attr_T>(p[f]);
		}
	}

	// /**
	//  * Get geogram facet attribute from ultimaille facet attribute
	//  */
	// template<typename MeshInput, int where, typename T>
	// constexpr void geo_attr_from_um_attr2(MeshInput &s, std::string attr_name, std::vector<T> &p, GEO::Mesh &m_out) {

	// 	using GEO_Attr_T = std::conditional_t<
	// 		std::is_same_v<T, double>, GEO::Numeric::float64,
	// 		std::conditional_t<std::is_same_v<T, int>, GEO::Numeric::uint32,
	// 		std::conditional_t<std::is_same_v<T, bool>, GEO::Numeric::uint8, T>>>;

   	// 	// auto geo_attr(
    //     // 	std::conditional_t<std::is_e<T, double>, decltype(m.vertices.attributes()), decltype(m.facets.attributes())>{}(), "flag"
	// 	// );

	// 	if constexpr(where == 1) {
	// 		GEO::Attribute<GEO_Attr_T> geo_attr(
	// 			m_out.facets.attributes(), attr_name
	// 		);

	// 		for (auto f : s.iter_facets()) {			
	// 			geo_attr[f] = static_cast<GEO_Attr_T>(p[f]);
	// 		}
	// 	}
	// }

	/**
	 * Get geogram facet attribute from ultimaille facet attribute
	 */
	// TODO refactor to avoid duplicate + change cell_facets, by anything
	template<typename T>
	void geo_attr_from_um_attr2(UM::Volume &s, std::string attr_name, std::vector<T> &p, GEO::Mesh &m) {

		using GEO_Attr_T = std::conditional_t<
			std::is_same_v<T, double>, GEO::Numeric::float64,
			std::conditional_t<std::is_same_v<T, int>, GEO::Numeric::uint32,
			std::conditional_t<std::is_same_v<T, bool>, GEO::Numeric::uint8, T>>>;

		GEO::Attribute<GEO_Attr_T> geo_attr(
			m.cell_facets.attributes(), attr_name
		);

		for (auto f : s.iter_facets()) {			
			geo_attr[f] = static_cast<GEO_Attr_T>(p[f]);
		}
	}

	template<int TWhere, typename TUMMesh, typename TAttr>
	void geo_attr_from_um_attr3(TUMMesh &s, std::string attr_name, std::vector<TAttr> &p, GEO::Mesh &m) {

		// TODO complete !
		using GEO_Attr_T = std::conditional_t<
			std::is_same_v<TAttr, double>, GEO::Numeric::float64,
			std::conditional_t<std::is_same_v<TAttr, int>, GEO::signed_index_t,
			std::conditional_t<std::is_same_v<TAttr, bool>, GEO::Numeric::uint8, TAttr>>>;

		// TODO complete !
		if constexpr (TWhere == GEO::MESH_CELL_FACETS) {
			GEO::Attribute<GEO_Attr_T> geo_attr(
				m.cell_facets.attributes(), attr_name
			);

			for (auto f : s.iter_facets()) {			
				geo_attr[f] = static_cast<GEO_Attr_T>(p[f]);
			}
		} 
		if constexpr (TWhere == GEO::MESH_FACETS) {
			GEO::Attribute<GEO_Attr_T> geo_attr(
				m.facets.attributes(), attr_name
			);

			for (auto f : s.iter_facets()) {			
				geo_attr[f] = static_cast<GEO_Attr_T>(p[f]);
			}
		}
	}

	template<int TWhere, typename TUMMesh, typename TAttr>
	void geo_attr_from_um_attr2(TUMMesh &m, std::shared_ptr<AttributeContainer<TAttr>> attr, std::string attr_name, GEO::Mesh &m_out) {
		geo_attr_from_um_attr3<TWhere>(m, attr_name, attr->data, m_out); 
	}

	/**
	 * Get UM mesh attribute from GEO mesh attribute
	 */
	template<int TWhere, typename TUMMesh, typename TAttr>
	void um_attr_from_geo_attr(GEO::Mesh &m, std::string attr_name, TUMMesh &m_out, std::shared_ptr<AttributeContainer<TAttr>> &um_attr) {

		// TODO complete !
		using GEO_Attr_T = std::conditional_t<
			std::is_same_v<TAttr, double>, GEO::Numeric::float64,
			std::conditional_t<std::is_same_v<TAttr, int>, GEO::signed_index_t,
			std::conditional_t<std::is_same_v<TAttr, bool>, GEO::Numeric::uint8, TAttr>>>;

		// TODO complete !
		if constexpr (TWhere == GEO::MESH_CELL_FACETS) {
			GEO::Attribute<GEO_Attr_T> geo_attr(
				m.cell_facets.attributes(), attr_name
			);

			for (auto f : m_out.iter_facets()) {		
				um_attr->data[f] = static_cast<TAttr>(geo_attr[f]);
			}
		} if constexpr (TWhere == GEO::MESH_FACETS) {
			GEO::Attribute<GEO_Attr_T> geo_attr(
				m.facets.attributes(), attr_name
			);

			for (auto f : m_out.iter_facets()) {		
				um_attr->data[f] = static_cast<TAttr>(geo_attr[f]);
			}
		}
	}

	// /**
	//  * Get ultimaille facet attribute from geogram facet attribute
	//  */
	// template<typename T>
	// void um_attr_from_geo_attr(GEO::Mesh &m, std::string attr_name, UM::Volume &s, UM::CellFacetAttribute<T> &um_attr) {

	// 	using GEO_Attr_T = std::conditional_t<
	// 		std::is_same_v<T, double>, GEO::Numeric::float64,
	// 		std::conditional_t<std::is_same_v<T, int>, GEO::Numeric::uint32,
	// 		std::conditional_t<std::is_same_v<T, bool>, GEO::Numeric::uint8, T>>>;

	// 	GEO::Attribute<GEO_Attr_T> geo_attr(
	// 		m.cell_facets.attributes(), attr_name
	// 	);

	// 	for (auto f : s.iter_facets()) {		
	// 		auto a = geo_attr[f];	
	// 		auto b = um_attr[f];	
	// 		um_attr[f] = static_cast<T>(geo_attr[f]);
	// 	}
	// }

	void combinatorial_update(GEO::Mesh &mesh, std::function<std::pair<SurfaceAttributes, VolumeAttributes>(TetBoundary&, SurfaceAttributes&, VolumeAttributes&)> f);
	
	// TODO add function to modify only attributes attributes_update
	// TODO add function to modify only geometry geometric_update



	// struct MeshBinding {

	// 	MeshBinding(GEO::Mesh &mesh);

	// 	GEO::Mesh &mesh;
	// 	UM::Tetrahedra tet;
	// 	TetBoundary tet_bound;

	// 	void update_tet();
	// 	void update_mesh();

	// 	template<typename TAttr>
	// 	void update_mesh_attr(std::shared_ptr<AttributeContainer<TAttr>> attr, std::string attr_name);
	// };

}