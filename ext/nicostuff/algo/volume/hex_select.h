#ifndef HEX_SELECT__H__
#define HEX_SELECT__H__

#include <ultimaille/all.h>
#include <framework/trace.h>


#include "toolbox.h"
#include <iostream>

#include "mesh_geom.h"


struct HexPointSelect :public PointAttribute<bool> {

	HexPointSelect(Hexahedra& hex, bool debug_lock = false) : hex(hex), PointAttribute<bool>(hex, false), debug_lock(debug_lock) {
		um_assert(hex.connected());
		debug_lock = false;
	}
	Hexahedra& hex;
	bool debug_lock;

	inline void save_state(PointAttribute<bool>& backup) {
		for (auto v : hex.iter_vertices()) backup[v] = (*this)[v];
	}
	inline void restore_state(PointAttribute<bool>& backup) {
		for (auto v : hex.iter_vertices()) (*this)[v] = backup[v];
	}

	

	inline void drop(std::string name = "select", bool inverse = false, bool forced = true) {
		if (debug_lock || forced) Drop(hex.points, static_cast<PointAttribute<bool>&> (*this))._skip_value(inverse).apply(name);
	}


	inline void fill(bool val) {
		for (auto v : hex.iter_vertices()) (*this)[v] = val;
		drop("fill",false,false);
	};

	inline void by_value(PointAttribute<int> attr, int attr_value, bool val) {
		for (auto v : hex.iter_vertices()) if (attr[v] == attr_value) (*this)[v] = val;
		drop("by value", false, false);
	};

	inline void by_facet_value(CellFacetAttribute<int>& attr, int attr_value, bool val) {
		for (auto h : hex.iter_halfedges())
			if (attr[h.facet()] == attr_value)
				(*this)[h.from()] = val;
		drop("by facet value", false, false);
	};

	inline void by_facet_non_value(CellFacetAttribute<int>& attr, int attr_value, bool val) {
		for (auto h : hex.iter_halfedges())
			if (attr[h.facet()] != attr_value)
				(*this)[h.from()] = val;
		drop("by_facet_non_value", false, false);
	};

	inline void invert() {
		for (auto v : hex.iter_vertices())
			(*this)[v] = !(*this)[v];
		drop("invert", false, false);
	}

	inline void dilate(bool val) {
		PointAttribute<bool> tmp(hex, !val);
		for (auto h : hex.iter_halfedges()) {
			if ((*this)[h.from()] == val) tmp[h.to()] = val;
		}
		for (auto v : hex.iter_vertices())
			if ((*this)[v] != val) (*this)[v] = tmp[v];
		drop("dilate", false, false);
	};




	inline void set_cell(int c, bool val) {
		FOR(lv, 8) (*this)[hex.vert(c, lv)] = val;
		drop("lock_cell", false, false);
	}





	inline void drop_cells(std::string name = "select", int min_nb_locked_vertices = 8) {
		Hexahedra select;
		ToolBox(select.points).copy_from(hex.points);
		for (auto c : hex.iter_cells()) {
			int nb_locked_vert = 0;
			FOR(lv, 8) if ((*this)[hex.vert(c, lv)])nb_locked_vert++;
			if (nb_locked_vert < min_nb_locked_vertices) continue;

			int nc = select.create_cells(1);
			FOR(lv, 8) select.vert(nc, lv) = hex.vert(c, lv);
		}
		select.delete_isolated_vertices();
		CellAttribute<double> SJ(select);
		ToolBox(select).eval_quality(SJ);
		Drop(select, SJ).apply(name);
	}

	inline void drop_cells_facet(CellFacetAttribute<int>& chart, std::string name = "select_facet", int min_nb_locked_vertices = 1) {
		Hexahedra select;
		CellFacetAttribute<int> chart_select(select);
		ToolBox(select.points).copy_from(hex.points);
		for (auto c : hex.iter_cells()) {
			int nb_locked_vert = 0;
			FOR(lv, 8) if ((*this)[hex.vert(c, lv)])nb_locked_vert++;
			if (nb_locked_vert < min_nb_locked_vertices) continue;
			int nc = select.create_cells(1);
			FOR(lv, 8) select.vert(nc, lv) = hex.vert(c, lv);
			FOR(lf, 6) chart_select[6 * nc + lf] = chart[6 * c + lf];
		}
		select.delete_isolated_vertices();
		Drop(select, chart_select)._shrink(0).apply(name);
	}

};
#endif