#pragma once

#include <ultimaille/all.h>

#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/GLUP/GLUP_private.h>
#include <geogram/basic/geometry.h>

#include "geom_ultimaille_binding.h"

namespace gl_draw {

	inline Quaternion align_with_uv(UM::vec3 u, UM::vec3 v) {
		v.normalize();
		u.normalize();
		
		if (u * v < -.99999) { 
			Quaternion res;  
			res.v = {1.0, 0., 0.}; 
			res.w = 0;
			return res; 
		}
		
		if (std::abs(u * v) > .99999)
			return Quaternion();

		UM::vec3 inbetwen_uv(v + u);
		inbetwen_uv.normalize();

		Quaternion res;
		res.w = v * inbetwen_uv; // scalar product with (1,0,0) divided by norm
		res.v = cross(inbetwen_uv, v); // cross product with (1,0,0) 
		
		return res;
	}

	static void draw_arrow(UM::vec3 A, UM::vec3 B, double r, int res, double cone_tube_ratio, GEO::vec4f color) {
		
		if (r == -1) 
			r = (B - A).norm() / 10.;

		UM::vec3 n = B - A;
		double l = std::max(1e-5, n.norm());
		n = n / l;

		Quaternion quat = align_with_uv({0, 0, 1}, n);
		auto M = quat.rotation_matrix();

		glupSetColor4fv(GLUP_FRONT_COLOR, color.data());
		auto mesh_width = glupGetMeshWidth();
		glupSetMeshWidth(0.);

		// Set temporary picking mode to constant to deactivate picking on GL primitive
		auto picking_mode = glupGetPickingMode();
		glupPickingMode(GLUP_PICK_CONSTANT);

		glupBegin(GLUP_TRIANGLES);

		for (int v = 0; v < res; v++) {
			// circle parametric equation
			auto pt = r * UM::vec3{cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0.};
			auto pt_nxt = r * UM::vec3{cos(2. * M_PI * double(v + 1) / double(res)), sin(2. * M_PI * double(v + 1) / double(res)), 0.};

			auto pt_base = (r * cone_tube_ratio) * UM::vec3{cos(2. * M_PI * double(v) / double(res)), sin(2. * M_PI * double(v) / double(res)), 0.};
			auto pt_nxt_base = (r * cone_tube_ratio) * UM::vec3{cos(2. * M_PI * double(v + 1) / double(res)), sin(2. * M_PI * double(v + 1) / double(res)), 0.};

			auto apex = UM::vec3{0, 0, l};

			auto h = UM::vec3{0, 0, std::max(0., l - 2.5 * r)};
			auto pt_cone_base = pt + h;
			auto pt_nxt_cone_base = pt_nxt + h;

			// Cone
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_cone_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_cone_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * apex + A).data());

			// Cone base
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_cone_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_cone_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * h + A).data());

			// Tube base
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * UM::vec3{0,0,0} + A).data());

			// Tube
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_nxt_base + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * (pt_nxt_base + h) + A).data());

			glupPrivateVertex3dv(um_bindings::geo_vec(M * (pt_nxt_base + h) + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * (pt_base + h) + A).data());
			glupPrivateVertex3dv(um_bindings::geo_vec(M * pt_base + A).data());

		}

		glupEnd();
		glupSetMeshWidth(mesh_width);
		glupPickingMode(picking_mode);
	}

	static void draw_path(std::vector<UM::vec3> &path, GEO::vec4f color, bool arrow) {

		auto picking_mode = glupGetPickingMode();
		glupPickingMode(GLUP_PICK_CONSTANT);

		glupSetColor4fv(GLUP_FRONT_COLOR, color.data());
		glupBegin(GLUP_POINTS);
		for (auto p : path) {
			double pp[3] = {p.x, p.y, p.z};
			glupPrivateVertex3dv(pp);
		}
		glupEnd();

		glupPickingMode(picking_mode);

		if (!arrow) {
			glupBegin(GLUP_LINES);
			for (int i = 0; i < path.size(); i+=2) {
				auto pA = path[i];
				auto pB = path[i + 1];
				double ppA[3] = {pA.x, pA.y, pA.z};
				double ppB[3] = {pB.x, pB.y, pB.z};
				glupPrivateVertex3dv(ppA);
				glupPrivateVertex3dv(ppB);
			}
			glupEnd();
		} else {
			for (int i = 0; i < path.size(); i+=2) {
				gl_draw::draw_arrow(path[i], path[i + 1], -1, 10, 0.5, color);
			}
		}

	}

	static void draw_path2(std::vector<UM::vec3> &path, GEO::vec4f color, bool arrow) {

		// auto picking_mode = glupGetPickingMode();
		// glupPickingMode(GLUP_PICK_CONSTANT);

		glupSetColor4fv(GLUP_FRONT_COLOR, color.data());

		glupBegin(GLUP_POINTS);
		for (auto p : path) {
			double pp[3] = {p.x, p.y, p.z};
			glupPrivateVertex3dv(pp);
		}
		glupEnd();

		// glupPickingMode(picking_mode);


		glupBegin(GLUP_LINES);
		for (int i = 0; i < static_cast<int>(path.size()) - 1; i++) {
			auto pA = path[i];
			auto pB = path[i + 1];
			double ppA[3] = {pA.x, pA.y, pA.z};
			double ppB[3] = {pB.x, pB.y, pB.z};
			glupPrivateVertex3dv(ppA);
			glupPrivateVertex3dv(ppB);
		}
		glupEnd();


	}

	static void draw_grid(double scale = 10., int n_lines = 9, double y = -.1) {

		glupSetColor4fv(GLUP_FRONT_COLOR, GEO::vec3f(0.2f,0.2f,0.2f).data());
		glupBegin(GLUP_LINES);
		for (int i = 0; i < n_lines; i++) {
			double j = i * (scale / n_lines) - scale * .5;
			double k = scale * .5;

			GEO::vec3 a(j, y, k);
			GEO::vec3 b(j, y, -k);
			
			glupPrivateVertex3dv(a.data());
			glupPrivateVertex3dv(b.data());
		}
		for (int i = 0; i < n_lines; i++) {
			double j = i * (scale / n_lines) - scale * .5;
			double k = scale * .5;

			GEO::vec3 a(k, y, j);
			GEO::vec3 b(-k, y, j);
			
			glupPrivateVertex3dv(a.data());
			glupPrivateVertex3dv(b.data());
		}
		glupEnd();
	}

	// TODO remove unit cell facet drawing !
	/**
	 * Draw overlay on cell facet
	 */
	static void draw_cell_facet_overlay(Mesh &mesh, index_t c, index_t lfacet, GEO::SimpleApplication::ColormapInfo colormap, GLUPdouble texCoord = 0., GLUPint thickness = 3.) {
		auto picking_mode = glupGetPickingMode();
		// glupPickingMode(GLUP_PICK_CONSTANT);
		// glupPickingId(-1);


		glupEnable(GLUP_TEXTURING);
		glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
		glBindTexture(GL_TEXTURE_2D, (GLuint) colormap.texture);
		glupTextureType(GLUP_TEXTURE_2D);
		glupTextureMode(GLUP_TEXTURE_REPLACE);
		glupPrivateTexCoord1d(texCoord);
		glupSetMeshWidth(thickness);

		glupBegin(GLUP_LINES);
		
		index_t facet_nb_verts = mesh.cells.facet_nb_vertices(c, lfacet);
		for(index_t lv = 0; lv < facet_nb_verts; lv++) {
			index_t v0_idx = mesh.cells.facet_vertex(c, lfacet, lv);
			index_t v1_idx = mesh.cells.facet_vertex(c, lfacet, (lv + 1) % facet_nb_verts);
			GEO::vec3 &a = mesh.vertices.point(v0_idx);
			GEO::vec3 &b = mesh.vertices.point(v1_idx);
			glupPrivateVertex3dv(a.data());
			glupPrivateVertex3dv(b.data());
		}

		glupDisable(GLUP_TEXTURING);
		glupEnd();
		
		// glupPickingId(idx);
		// glupPickingMode(picking_mode);

		// glupBegin(GLUP_POINTS);

		// for(index_t lv = 0; lv < facet_nb_verts; lv++) {
		// 	index_t v0_idx = mesh.cells.facet_vertex(c, lfacet, lv);
		// 	GEO::vec3 &a = mesh.vertices.point(v0_idx);
		// 	glupPrivateVertex3dv(a.data());
		// }


		// glupEnd();
	}

	// TODO remove unit cell drawing !
	/**
	 * Draw overlay on cell
	 */
	static void draw_cell_overlay(Mesh &mesh, index_t c, GEO::SimpleApplication::ColormapInfo colormap, GLUPdouble texCoord = 0., GLUPint thickness = 3.) {
		auto picking_mode = glupGetPickingMode();
		glupPickingMode(GLUP_PICK_CONSTANT);
		int idx = glupGetPickingId();
		glupPickingId(idx);

		glupEnable(GLUP_TEXTURING);
		glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
		glBindTexture(GL_TEXTURE_2D, (GLuint) colormap.texture);
		glupTextureType(GLUP_TEXTURE_2D);
		glupTextureMode(GLUP_TEXTURE_REPLACE);
		glupPrivateTexCoord1d(texCoord);
		glupSetMeshWidth(thickness);

		glupBegin(GLUP_LINES);

		for(index_t e = 0; e < mesh.cells.nb_edges(c); e++) {
			index_t v0_idx = mesh.cells.edge_vertex(c, e, 0);
			index_t v1_idx = mesh.cells.edge_vertex(c, e, 1);
			GEO::vec3 &a = mesh.vertices.point(v0_idx);
			GEO::vec3 &b = mesh.vertices.point(v1_idx);
			glupPrivateVertex3dv(a.data());
			glupPrivateVertex3dv(b.data());
		}
		
		glupDisable(GLUP_TEXTURING);
		glupEnd();

		glupPickingId(idx);
		glupPickingMode(picking_mode);


	}

	/**
	 * Draw overlay on cells
	 */
	static void draw_cells_overlay(Mesh &mesh, std::vector<int> cells, GEO::SimpleApplication::ColormapInfo colormap, GLUPdouble texCoord = 0., GLUPint thickness = 3.) {
		auto picking_mode = glupGetPickingMode();
		glupPickingMode(GLUP_PICK_CONSTANT);
		int idx = glupGetPickingId();
		glupPickingId(idx);

		glupEnable(GLUP_TEXTURING);
		glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
		glBindTexture(GL_TEXTURE_2D, (GLuint) colormap.texture);
		glupTextureType(GLUP_TEXTURE_2D);
		glupTextureMode(GLUP_TEXTURE_REPLACE);
		glupPrivateTexCoord1d(texCoord);
		glupSetMeshWidth(thickness);

		glupBegin(GLUP_LINES);

		for (auto c : cells) {
			for(index_t e = 0; e < mesh.cells.nb_edges(c); e++) {
				index_t v0_idx = mesh.cells.edge_vertex(c, e, 0);
				index_t v1_idx = mesh.cells.edge_vertex(c, e, 1);
				GEO::vec3 &a = mesh.vertices.point(v0_idx);
				GEO::vec3 &b = mesh.vertices.point(v1_idx);
				glupPrivateVertex3dv(a.data());
				glupPrivateVertex3dv(b.data());
			}
		}
		glupDisable(GLUP_TEXTURING);
		glupEnd();

		glupPickingId(idx);
		glupPickingMode(picking_mode);


	}

	static void draw_cell_facets_overlay(Mesh &mesh, std::vector<std::pair<int, int>> cell_lfacets, GEO::SimpleApplication::ColormapInfo colormap, GLUPdouble texCoord = 0., GLUPint thickness = 3.) {
		auto picking_mode = glupGetPickingMode();
		// glupPickingMode(GLUP_PICK_CONSTANT);
		// glupPickingId(-1);


		glupEnable(GLUP_TEXTURING);
		glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
		glBindTexture(GL_TEXTURE_2D, (GLuint) colormap.texture);
		glupTextureType(GLUP_TEXTURE_2D);
		glupTextureMode(GLUP_TEXTURE_REPLACE);
		glupPrivateTexCoord1d(texCoord);
		glupSetMeshWidth(thickness);

		glupBegin(GLUP_LINES);
		
		for (auto [c, lfacet] : cell_lfacets) {
			index_t facet_nb_verts = mesh.cells.facet_nb_vertices(c, lfacet);
			for(index_t lv = 0; lv < facet_nb_verts; lv++) {
				index_t v0_idx = mesh.cells.facet_vertex(c, lfacet, lv);
				index_t v1_idx = mesh.cells.facet_vertex(c, lfacet, (lv + 1) % facet_nb_verts);
				GEO::vec3 &a = mesh.vertices.point(v0_idx);
				GEO::vec3 &b = mesh.vertices.point(v1_idx);
				glupPrivateVertex3dv(a.data());
				glupPrivateVertex3dv(b.data());
			}
		}

		glupDisable(GLUP_TEXTURING);
		glupEnd();
	}

	static void test2(Mesh &mesh, std::vector<int> cells, GEO::SimpleApplication::ColormapInfo colormap, GLUPdouble texCoord = 0., GLUPint thickness = 3.) {

		if (cells.size() <= 0)
			return;

		index_t fc = *mesh.cells.begin();
		index_t facet_nb_verts = mesh.cells.facet_nb_vertices(fc, 0);

		GLUPprimitive p = GLUP_QUADS;
		switch (facet_nb_verts)
		{
		case 3:
			p = GLUP_TRIANGLES;
			break;
		case 4:
			p = GLUP_QUADS;
			break;
		default:
			throw new std::runtime_error("Weird ! facet has more of 4 vertices or less of 3...");
			break;
		}

		// auto picking_mode = glupGetPickingMode();
		// glupPickingMode(GLUP_PICK_CONSTANT);
		// glupPickingId(-1);

		glupEnable(GLUP_TEXTURING);
		glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
		glBindTexture(GL_TEXTURE_2D, (GLuint) colormap.texture);
		glupTextureType(GLUP_TEXTURE_2D);
		glupTextureMode(GLUP_TEXTURE_REPLACE);
		glupPrivateTexCoord1d(texCoord);
		glupSetMeshWidth(thickness);

		glupBegin(p);

		for (auto c : cells) {
			for (int lf = 0; lf < mesh.cells.nb_facets(c); lf++) {
				index_t f = mesh.cells.facet(c, lf);
				if (mesh.cell_facets.adjacent_cell(f) != NO_CELL) {
					for(index_t lv = 0; lv < facet_nb_verts; lv++) {
						index_t v = mesh.cells.facet_vertex(c, lf, lv);
						GEO::vec3 &p = mesh.vertices.point(v);
						glupPrivateVertex3dv(p.data());
					}
				}
			}
		}


		glupDisable(GLUP_TEXTURING);
		glupEnd();
	}

	static void test(Mesh &mesh, std::vector<std::pair<int, int>> cell_lfacets, GEO::SimpleApplication::ColormapInfo colormap, GLUPdouble texCoord = 0., GLUPint thickness = 3.) {

		if (cell_lfacets.size() <= 0)
			return;

		// auto picking_mode = glupGetPickingMode();
		// glupPickingMode(GLUP_PICK_CONSTANT);
		// glupPickingId(-1);

		glupEnable(GLUP_TEXTURING);
		glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
		glBindTexture(GL_TEXTURE_2D, (GLuint) colormap.texture);
		glupTextureType(GLUP_TEXTURE_2D);
		glupTextureMode(GLUP_TEXTURE_REPLACE);
		glupPrivateTexCoord1d(texCoord);
		glupSetMeshWidth(thickness);


		auto [c, lfacet] = *cell_lfacets.begin();

		index_t facet_nb_verts = mesh.cells.facet_nb_vertices(c, lfacet);

		GLUPprimitive p = GLUP_QUADS;
		switch (facet_nb_verts)
		{
		case 3:
			p = GLUP_TRIANGLES;
			break;
		case 4:
			p = GLUP_QUADS;
			break;
		default:
			throw new std::runtime_error("Weird ! facet has more of 4 vertices or less of 3...");
			break;
		}


		glupBegin(p);
		
		for (auto [c, lfacet] : cell_lfacets) {
			for(index_t lv = 0; lv < facet_nb_verts; lv++) {
				index_t v = mesh.cells.facet_vertex(c, lfacet, lv);
				GEO::vec3 &p = mesh.vertices.point(v);
				glupPrivateVertex3dv(p.data());
			}
		}

		glupDisable(GLUP_TEXTURING);
		glupEnd();
	}

	static void draw_axis(GEO::vec3f origin = GEO::vec3f(0, 0, 0)) {
		glupSetColor4fv(GLUP_FRONT_COLOR, (origin + GEO::vec3f(1.f,0.f,0.f)).data());
		glupBegin(GLUP_LINES);
		glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
		glupPrivateVertex3dv(GEO::vec3(1, 0, 0).data());
		glupEnd();

		glupSetColor4fv(GLUP_FRONT_COLOR, (origin + GEO::vec3f(0.f,1.f,0.f)).data());
		glupBegin(GLUP_LINES);
		glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
		glupPrivateVertex3dv(GEO::vec3(0, 1, 0).data());
		glupEnd();

		glupSetColor4fv(GLUP_FRONT_COLOR, (origin + GEO::vec3f(0.f,0.f,1.f)).data());
		glupBegin(GLUP_LINES);
		glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
		glupPrivateVertex3dv(GEO::vec3(0, 0, 1).data());
		glupEnd();

		glupSetColor4fv(GLUP_FRONT_COLOR, (origin + GEO::vec3f(1.f, 1.f, 1.f)).data());
		glupBegin(GLUP_POINTS);
		glupPrivateVertex3dv(GEO::vec3(0, 0, 0).data());
		glupEnd();
	}

}