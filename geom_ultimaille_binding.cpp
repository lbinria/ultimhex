#include <functional>

#include <ultimaille/all.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/geofile.h>

#include "geom_ultimaille_binding.h"

namespace um_bindings {

	template<typename T>
	void set_um_vertices_from_geo_vertices(GEO::Mesh &m, T &m_out) {
		for (auto v : m.vertices) {
			auto p = m.vertices.point(v);
			m_out.points[v] = UM::vec3{p.x, p.y, p.z};
		}
	}

	template<typename T>
	void set_geo_vertices_from_um_vertices(T &m, GEO::Mesh &m_out) {
		for (auto v : m.iter_vertices()) {
			m_out.vertices.point(v).pos() = v.pos();
		}
	}

	void um_tri_from_geo_mesh(GEO::Mesh &m, UM::Triangles &m_out) {
		m_out.clear();

		// TODO here not true ! m.vertices.nb ! should keep only surface points
		m_out.points.create_points(m.vertices.nb());
		m_out.create_facets(m.facets.nb());

		set_um_vertices_from_geo_vertices(m, m_out);

		for (auto f : m.facets) {
			m_out.vert(f, 0) = m.facets.vertex(f, 0);
			m_out.vert(f, 1) = m.facets.vertex(f, 1);
			m_out.vert(f, 2) = m.facets.vertex(f, 2);
		}
	}

	void um_tet_from_geo_mesh(GEO::Mesh &m, UM::Tetrahedra &m_out) {
		m_out.clear();

		m_out.points.create_points(m.vertices.nb());
		m_out.create_cells(m.cells.nb());

		set_um_vertices_from_geo_vertices(m, m_out);

		for (auto c : m.cells) {
			m_out.vert(c, 0) = m.cells.vertex(c, 0);
			m_out.vert(c, 1) = m.cells.vertex(c, 1);
			m_out.vert(c, 2) = m.cells.vertex(c, 2);
			m_out.vert(c, 3) = m.cells.vertex(c, 3);
		}
	}

	void um_hex_from_geo_mesh(GEO::Mesh &m, UM::Hexahedra &m_out) {
		m_out.clear();

		m_out.points.create_points(m.vertices.nb());
		m_out.create_cells(m.cells.nb());

		set_um_vertices_from_geo_vertices(m, m_out);

		for (auto c : m.cells) {
			m_out.vert(c, 0) = m.cells.vertex(c, 0);
			m_out.vert(c, 1) = m.cells.vertex(c, 1);
			m_out.vert(c, 2) = m.cells.vertex(c, 2);
			m_out.vert(c, 3) = m.cells.vertex(c, 3);
			m_out.vert(c, 4) = m.cells.vertex(c, 4);
			m_out.vert(c, 5) = m.cells.vertex(c, 5);
			m_out.vert(c, 6) = m.cells.vertex(c, 6);
			m_out.vert(c, 7) = m.cells.vertex(c, 7);
		}
	}

	// Check !
	void geo_mesh_from_um_tri(UM::Triangles &m, GEO::Mesh &m_out, bool clear) {

		if (clear)
			m_out.clear(false, false);

		int v_off = m_out.vertices.create_vertices(m.nverts());
		int f_off = m_out.facets.create_facets(m.nfacets(), 3);

		for (auto v : m.iter_vertices()) {
			auto p = v.pos();
			m_out.vertices.point(v_off + v) = GEO::vec3(p.x, p.y, p.z);
		}

		for (auto f : m.iter_facets()) {
			m_out.facets.set_vertex(f + f_off, 0, m.vert(f, 0) + v_off);
			m_out.facets.set_vertex(f + f_off, 1, m.vert(f, 1) + v_off);
			m_out.facets.set_vertex(f + f_off, 2, m.vert(f, 2) + v_off);
		}

	}

	void geo_mesh_from_um_tet(UM::Tetrahedra &m, GEO::Mesh &m_out, bool clear) {
		assert(m.connected());

		if (clear)
			m_out.clear(false, false);

		int v_off = m_out.vertices.create_vertices(m.nverts());
		int c_off = m_out.cells.create_cells(m.ncells(), GEO::MESH_TET);

		for (auto v : m.iter_vertices()) {
			auto p = v.pos();
			m_out.vertices.point(v_off + v) = GEO::vec3(p.x, p.y, p.z);
		}		

		for (auto c : m.iter_cells()) {
			m_out.cells.set_vertex(c + c_off, 0, m.vert(c, 0) + v_off);
			m_out.cells.set_vertex(c + c_off, 1, m.vert(c, 1) + v_off);
			m_out.cells.set_vertex(c + c_off, 2, m.vert(c, 2) + v_off);
			m_out.cells.set_vertex(c + c_off, 3, m.vert(c, 3) + v_off);
		}
		
	}

	void geo_mesh_from_tetboundary(TetBoundary &m, GEO::Mesh &m_out) {

		m_out.clear();

		int v_off = m_out.vertices.create_vertices(m.tet.nverts());
		int f_off = m_out.facets.create_facets(m.tri.nfacets(), 3);
		int c_off = m_out.cells.create_cells(m.tet.ncells(), GEO::MESH_TET);

		for (auto v : m.tet.iter_vertices()) {
			auto p = v.pos();
			m_out.vertices.point(v_off + v) = GEO::vec3(p.x, p.y, p.z);
		}

		for (auto f : m.tri.iter_facets()) {
			m_out.facets.set_vertex(f + f_off, 0, m.tri.vert(f, 0) + v_off);
			m_out.facets.set_vertex(f + f_off, 1, m.tri.vert(f, 1) + v_off);
			m_out.facets.set_vertex(f + f_off, 2, m.tri.vert(f, 2) + v_off);
		}

		for (auto c : m.tet.iter_cells()) {
			m_out.cells.set_vertex(c + c_off, 0, m.tet.vert(c, 0) + v_off);
			m_out.cells.set_vertex(c + c_off, 1, m.tet.vert(c, 1) + v_off);
			m_out.cells.set_vertex(c + c_off, 2, m.tet.vert(c, 2) + v_off);
			m_out.cells.set_vertex(c + c_off, 3, m.tet.vert(c, 3) + v_off);
		}
	}

	void geo_mesh_from_um_quad(UM::Quads &m, GEO::Mesh &m_out, bool clear) {

		if (clear)
			m_out.clear(false, false);

		int v_off = m_out.vertices.create_vertices(m.nverts());
		int f_off = m_out.facets.create_facets(m.nfacets(), 4);

		for (auto v : m.iter_vertices()) {
			auto p = v.pos();
			m_out.vertices.point(v_off + v) = GEO::vec3(p.x, p.y, p.z);
		}

		for (auto f : m.iter_facets()) {
			m_out.facets.set_vertex(f + f_off, 0, m.vert(f, 0) + v_off);
			m_out.facets.set_vertex(f + f_off, 1, m.vert(f, 1) + v_off);
			m_out.facets.set_vertex(f + f_off, 2, m.vert(f, 2) + v_off);
			m_out.facets.set_vertex(f + f_off, 3, m.vert(f, 3) + v_off);
		}

	}

	void geo_mesh_from_hexboundary(HexBoundary &m, GEO::Mesh &m_out) {

		m_out.clear(true, true);

		int v_off = m_out.vertices.create_vertices(m.hex.nverts());
		int f_off = m_out.facets.create_facets(m.quad.nfacets(), 4);
		int c_off = m_out.cells.create_cells(m.hex.ncells(), GEO::MESH_HEX);

		for (auto v : m.hex.iter_vertices()) {
			auto p = v.pos();
			m_out.vertices.point(v_off + v) = GEO::vec3(p.x, p.y, p.z);
		}

		for (auto f : m.quad.iter_facets()) {
			m_out.facets.set_vertex(f + f_off, 0, m.quad.vert(f, 0) + v_off);
			m_out.facets.set_vertex(f + f_off, 1, m.quad.vert(f, 1) + v_off);
			m_out.facets.set_vertex(f + f_off, 2, m.quad.vert(f, 2) + v_off);
			m_out.facets.set_vertex(f + f_off, 3, m.quad.vert(f, 3) + v_off);
		}

		for (auto c : m.hex.iter_cells()) {
			m_out.cells.set_vertex(c + c_off, 0, m.hex.vert(c, 0) + v_off);
			m_out.cells.set_vertex(c + c_off, 1, m.hex.vert(c, 1) + v_off);
			m_out.cells.set_vertex(c + c_off, 2, m.hex.vert(c, 2) + v_off);
			m_out.cells.set_vertex(c + c_off, 3, m.hex.vert(c, 3) + v_off);
			m_out.cells.set_vertex(c + c_off, 4, m.hex.vert(c, 4) + v_off);
			m_out.cells.set_vertex(c + c_off, 5, m.hex.vert(c, 5) + v_off);
			m_out.cells.set_vertex(c + c_off, 6, m.hex.vert(c, 6) + v_off);
			m_out.cells.set_vertex(c + c_off, 7, m.hex.vert(c, 7) + v_off);
		}
	}

	void geo_mesh_from_um_hex(UM::Hexahedra &m, GEO::Mesh &m_out, bool clear) {
		// TODO see if necessary ?
		assert(m.connected());

		if (clear)
			m_out.clear(false, false);

		int v_off = m_out.vertices.create_vertices(m.nverts());
		int c_off = m_out.cells.create_cells(m.ncells(), GEO::MESH_HEX);

		for (auto v : m.iter_vertices()) {
			auto p = v.pos();
			m_out.vertices.point(v_off + v) = GEO::vec3(p.x, p.y, p.z);
		}

		for (auto c : m.iter_cells()) {
			m_out.cells.set_vertex(c + c_off, 0, m.vert(c, 0) + v_off);
			m_out.cells.set_vertex(c + c_off, 1, m.vert(c, 1) + v_off);
			m_out.cells.set_vertex(c + c_off, 2, m.vert(c, 2) + v_off);
			m_out.cells.set_vertex(c + c_off, 3, m.vert(c, 3) + v_off);
			m_out.cells.set_vertex(c + c_off, 4, m.vert(c, 4) + v_off);
			m_out.cells.set_vertex(c + c_off, 5, m.vert(c, 5) + v_off);
			m_out.cells.set_vertex(c + c_off, 6, m.vert(c, 6) + v_off);
			m_out.cells.set_vertex(c + c_off, 7, m.vert(c, 7) + v_off);
		}
		
	}

	std::vector<std::string> split_string(const std::string& str, char delimiter) {
		std::vector<std::string> result;
		std::stringstream ss(str);
		std::string token;

		while (std::getline(ss, token, delimiter)) {
			result.push_back(token);
		}

		return result;
	}

	// void combinatorial_update(GEO::Mesh &mesh, std::function<std::pair<SurfaceAttributes, VolumeAttributes>(TetBoundary&, SurfaceAttributes&, VolumeAttributes&)> f) {
	// 	Tetrahedra tet;
	// 	um_tet_from_geo_mesh(mesh, tet);

	// 	// TODO get attribute from mesh to SurfaceAttributes and VolumeAttributes
	// 	SurfaceAttributes input_surf_attr = {};
	// 	VolumeAttributes input_vol_attr = {};

	// 	GEO::vector<std::string> cell_facets_attr_names;
	// 	mesh.cell_facets.attributes().list_attribute_names(cell_facets_attr_names);


   	// 	std::vector<std::string> parts = split_string(mesh.get_attributes(), ';');
	// 	for (auto p : parts) {
	// 		GEO::MeshElementsFlags where;
	// 		std::string attr_name;
	// 		GEO::index_t component;
	// 		mesh.parse_attribute_name(p, where, attr_name, component);

	// 		if (where == GEO::MeshElementsFlags::MESH_VERTICES) {
				
	// 			GEO::AttributeStore *store = mesh.vertices.attributes().find_attribute_store(attr_name);
	// 			if (store->elements_type_matches(typeid(GEO::vec3).name())) {
	// 				// GEO::Attribute<GEO::vec3> geo_attr(
	// 				// 	mesh.vertices.attributes(), attr_name
	// 				// );
	// 				// std::cout << "ok" << std::endl;

	// 			}




	// 		} if (where == GEO::MeshElementsFlags::MESH_CELL_FACETS) {

	// 			GEO::AttributeStore *store = mesh.cell_facets.attributes().find_attribute_store(attr_name);
	// 			if (store->elements_type_matches(typeid(GEO::Numeric::uint32).name())) {
	// 				CellFacetAttribute<int> um_attr(tet);
	// 				um_attr_from_geo_attr(mesh, attr_name, tet, um_attr);
	// 				input_vol_attr.cell_facets.push_back({attr_name, um_attr.ptr});
	// 			}

	// 		}
	// 	}

	// 	tet.connect();
	// 	TetBoundary tet_bound(tet);
	// 	auto [output_surf_attr, output_vol_attr] = f(tet_bound, input_surf_attr, input_vol_attr);
	// 	geo_mesh_from_um_tet(tet_bound.tet, mesh);
	// 	geo_mesh_from_um_tri(tet_bound.tri, mesh, false);

	// 	// Transfert attributes from UM mesh to GEO mesh

	// 	// TODO treat other than facets ! and cells facet in volume
	// 	// Attributes
	// 	std::vector<NamedContainer> A[3] = {output_surf_attr.points, output_surf_attr.facets, output_surf_attr.corners};

	// 	for (int z=0; z<3; z++) {
	// 		auto &att = A[z];

	// 		for (int i=0; i<static_cast<int>(att.size()); i++) {
	// 			std::string name = att[i].first;
	// 			std::shared_ptr<GenericAttributeContainer> ptr = att[i].second;

	// 			// TODO encapsulate in function
	// 			if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tri, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<double>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tri, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tri, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec3>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tri, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<bool>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tri, name, cont_ptr->data, mesh);
	// 			} else {
	// 				assert(false);
	// 			}
	// 		}
	// 	}

	// 	// TODO refactor to avoid duplicate code
	// 	// Attributes
	// 	std::vector<NamedContainer> VA[4] = {output_vol_attr.points, output_vol_attr.cells, output_vol_attr.cell_facets, output_vol_attr.cell_corners};

	// 	for (int z=0; z<4; z++) {
	// 		auto &att = VA[z];

	// 		for (int i=0; i<static_cast<int>(att.size()); i++) {
	// 			std::string name = att[i].first;
	// 			std::shared_ptr<GenericAttributeContainer> ptr = att[i].second;

	// 			if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tet, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<double>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tet, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tet, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec3>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tet, name, cont_ptr->data, mesh);
	// 			} else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<bool>>(ptr); cont_ptr.get()!=nullptr) {
	// 				geo_attr_from_um_attr2(tet_bound.tet, name, cont_ptr->data, mesh);
	// 			} else {
	// 				assert(false);
	// 			}
	// 		}
	// 	}
	// }


	// MeshBinding::MeshBinding(GEO::Mesh &mesh) : mesh(mesh), tet(Tetrahedra()), tet_bound(tet) {

	// }

	// void MeshBinding::update_tet() {
		
	// 	// Disconnect !
	// 	if (tet_bound.tet.connected())
	// 		tet_bound.tet.disconnect();

	// 	// Recreate tet from scratch
	// 	um_tet_from_geo_mesh(mesh, tet_bound.tet);

	// 	// Reconnect for update
	// 	// TODO see if possible to do this without connecting
	// 	tet_bound.tet.connect();

	// 	// Update tri of tet boundaries
	// 	tet_bound.update();
	// }

	// void MeshBinding::update_mesh() {
	// 	geo_mesh_from_um_tet(tet_bound.tet, mesh);
	// }

	// template<typename TAttr>
	// void MeshBinding::update_mesh_attr(std::shared_ptr<AttributeContainer<TAttr>> attr, std::string attr_name) {
	// 	geo_attr_from_um_attr2(tet_bound.tet, attr_name, attr->data, mesh); 
	// }

}