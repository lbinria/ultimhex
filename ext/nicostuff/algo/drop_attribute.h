#ifndef DROP__ATTRIBUTE__H__
#define DROP__ATTRIBUTE__H__


#include <ultimaille/all.h>
#include <basic.h>

#include <algo/toolbox.h>
#include <algo/drop_mesh.h>
#include <algo/glyphs.h>


void extract_iso(Tetrahedra& m, Triangles& tr, FacetAttribute<double>& attr, CellCornerAttribute<double>& scal, double iso);

struct DropPolyLineGeometry {
    DropPolyLineGeometry(PolyLine& m);
    PolyLine& m;
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    DropPolyLineGeometry& _force_radius(double radius, bool relative_to_bbox = false) ;
    double __forced_radius;

    void auto_set_radius();

    void apply_arrow(std::string name = "");

    void apply_capsule(std::string name = "");

    void apply_disk(std::string name = "");
};



//    _____        _           _       _____        _   
//   |  __ \      (_)         | |     / ____|      | |  
//   | |__) |__    _   _ __   | |_   | (___   ___  | |_ 
//   |  ___/ _ \  | | | '_ \  | __|   \___ \ / _ \ | __|
//   | |  | (_) | | | | | | | | |_    ____) |  __/ | |_ 
//   |_|   \___/  |_| |_| |_|  \__|  |_____/ \___|  \__|
//                                                            

template<class T>
struct Drop<PointSet, PointAttribute<T>> {
    Drop(PointSet& pts, PointAttribute<T>& attr) : pts(pts), attr(attr) {}
    PointSet& pts;
    PointAttribute<T>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(T, skip_value, T(-1));
    ARG(std::string, just_save_filename, "");

    void apply(std::string name = "") {
        PointSet ptsout;
        PointAttribute<T> ptsout_attr(ptsout);

        FOR(v, pts.size()) if (attr[v] != __skip_value && __select(v)) {
            int nv = ptsout.create_points(1);
            ptsout_attr[nv] = attr[v];
            ptsout[nv] = pts[v];
        }
        DropPointSet(ptsout).add(ptsout_attr, "attr")._just_save_filename(__just_save_filename)._active_point_attribute("attr").apply(name);
    }
};


template<> struct Drop<PointSet, PointAttribute<vec3> > {
    Drop(PointSet& pts, PointAttribute<vec3>& attr) : pts(pts), attr(attr), __forced_length(0), __forced_radius(-1) {}
    PointSet& pts;
    PointAttribute<vec3>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_length(double length, bool relative_to_bbox = false) {
        __forced_length = length;
        if (relative_to_bbox) {
            auto bbox = ToolBox(pts).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_length *= size;
        }
        return *this;
    };
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(pts).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    double __forced_length;
    void apply_arrow(std::string name = "") {
        apply(name, false);
    }
    void apply_disk(std::string name = "") {
        apply(name, true);
    }
    void apply(std::string name = "", bool use_disk = false, bool use_line = true) {
        PolyLine pl;
        int nb_pts = 0;
        FOR(v, pts.size()) if (__select(v)) nb_pts++;
        pl.create_edges(nb_pts);
        pl.points.create_points(2 * nb_pts);
        int e = 0;
        FOR(v, pts.size()) if (__select(v)) {
            vec3 vect = attr[v];
            if (__forced_length != 0 && vect.norm() > 1e-20)
                vect = __forced_length * vect.normalized();
            FOR(i, 2) pl.vert(e, i) = 2 * e + i;
            pl.points[2 * e] = pts[v];
            pl.points[2 * e + 1] = pts[v] + vect;
            e++;
        }

        if (use_disk) DropPolyLineGeometry(pl)._just_save_filename(__just_save_filename)._resolution(__resolution)._force_radius(__forced_radius).apply_disk(name);
        else if (use_line) DropPolyLine(pl)._just_save_filename(__just_save_filename).apply(name);
        else DropPolyLineGeometry(pl)._just_save_filename(__just_save_filename)._resolution(__resolution)._force_radius(__forced_radius).apply_arrow(name);


    }
};


//    _____        _               _        _                
//   |  __ \      | |             | |      (_)               
//   | |__) |__   | |  _   _      | |       _   _ __     ___ 
//   |  ___/ _ \  | | | | | |     | |      | | | '_ \   / _ \
//   | |  | (_) | | | | |_| |     | |____  | | | | | | |  __/
//   |_|   \___/  |_|  \__, |     |______| |_| |_| |_|  \___|
//                      __/ |                                
//                     |___/                                 

template <typename T> concept PolyLineConceptBase = std::is_base_of<PolyLine, T>::value;


template<PolyLineConceptBase C,class T>
struct Drop<C, PointAttribute<T>> {
    Drop(C& m, PointAttribute<T>& attr) : m(m), attr(attr), __forced_radius(-1) {}
    C& m;
    PointAttribute<T>& attr;
    ARG(T, skip_value, T(-1));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    void auto_set_radius() {
        if (__forced_radius != -1) return;
        vec2 ave_edge_length(0, 0);
        FOR(s, m.nedges())
            if ((m.points[m.vert(s, 0)] - m.points[m.vert(s, 1)]).norm2() > 1e-10)
                ave_edge_length += vec2((m.points[m.vert(s, 0)] - m.points[m.vert(s, 1)]).norm(), 1);
        __forced_radius = .01 * ave_edge_length[0] / ave_edge_length[1];
    }
    void apply_arrow(std::string name = "") {
        auto_set_radius();
        Polygons outm;
        PointAttribute<T> out_attr(outm);
        FOR(s, m.nedges())
            if (__select(s) && attr[m.vert(s, 0)] != __skip_value && attr[m.vert(s, 1)] != __skip_value)
                Glyphs::Arrow(outm).resolution(__resolution).radius(__forced_radius).apply(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)], out_attr, attr[m.vert(s, 0)], attr[m.vert(s, 1)]);
        DropSurface(outm).add(out_attr, "out_attr")
            ._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false)._active_point_attribute("out_attr").apply(name);
    }
    void apply_capsule(std::string name = "") {
        auto_set_radius();
        Polygons outm;
        PointAttribute<T> out_attr(outm);
        FOR(s, m.nedges())
            if (__select(s) && attr[m.vert(s, 0)] != __skip_value && attr[m.vert(s, 1)] != __skip_value)
                Glyphs::Capsule(outm).resolution(__resolution).radius(__forced_radius).apply(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)], out_attr, attr[m.vert(s, 0)], attr[m.vert(s, 1)]);
        DropSurface(outm).add(out_attr, "out_attr")
            ._just_save_filename(__just_save_filename)._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false)._active_point_attribute("out_attr").apply(name);
    }
    void apply_wireframe(std::string name = "") {
        PolyLine pl;
        ToolBox(pl.points).copy_from(m.points);
        
        PointAttribute<T> out_attr(pl);
        FOR(s, m.nedges()) {
            T val = attr[s];
            if (!__select(s) || val == __skip_value) continue;
            int spl = pl.create_edges(1);
            FOR(lv, 2) pl.vert(spl, lv) = m.vert(s, lv);
            FOR(lv, 2) out_attr[pl.vert(spl, lv)] = attr[m.vert(s, lv)];
        }

        DropPolyLine(pl)._just_save_filename(__just_save_filename).add(out_attr, "out_attr")._active_point_attribute("out_attr").apply(name);
    }
};



template<PolyLineConceptBase C, class T>
struct Drop<C, EdgeAttribute<T>> {
    Drop(C& m, EdgeAttribute<T>& attr) : m(m), attr(attr), __forced_radius(-1) {}
    C& m;
    EdgeAttribute<T>& attr;
    ARG(T, skip_value, T(-1));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    void auto_set_radius() {
        if (__forced_radius != -1) return;
        vec2 ave_edge_length(0, 0);
        FOR(s, m.nedges())
            if ((m.points[m.vert(s, 0)] - m.points[m.vert(s, 1)]).norm2() > 1e-10)
                ave_edge_length += vec2((m.points[m.vert(s, 0)] - m.points[m.vert(s, 1)]).norm(), 1);
        __forced_radius = .01 * ave_edge_length[0] / ave_edge_length[1];
    }

    ARG(int, render_style, 0); // see code below 
    void apply(std::string name = "") {
        switch (__render_style) {
        case 0: apply_wireframe(name); break;
        case 1: apply_arrow(name); break;
        case 2: apply_capsule(name); break;
        }
    }


    void apply_arrow(std::string name = "") {
        auto_set_radius();
        Polygons outm;
        FacetAttribute<T> out_attr(outm);
        FOR(s, m.nedges()) {
            T attr_s = attr[s];

            if (__select(s) && attr[s] != __skip_value)
                Glyphs::Arrow(outm).resolution(__resolution).radius(__forced_radius)
                .apply(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)], out_attr, attr_s);
        }
        DropSurface(outm)._just_save_filename(__just_save_filename).add(out_attr, "out_attr")
            ._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false)._active_facet_attribute("out_attr").apply(name);
    }
    void apply_capsule(std::string name = "") {
        auto_set_radius();
        Polygons outm;
        FacetAttribute<T> out_attr(outm);
        FOR(s, m.nedges()) {
            T attr_s = attr[s];
            if (__select(s) && attr_s != __skip_value)
                Glyphs::Capsule(outm).resolution(__resolution)
                .radius(__forced_radius).apply<T>(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)], out_attr, attr_s);
        }
        DropSurface(outm)._just_save_filename(__just_save_filename).add(out_attr, "out_attr")
            ._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false)._active_facet_attribute("out_attr").apply(name);
    }
    void apply_wireframe(std::string name = "") {
        PolyLine pl;
        ToolBox(pl.points).copy_from(m.points);

        EdgeAttribute<T> out_attr(pl);
        FOR(s, m.nedges()) {
            T val = attr[s];
            if (!__select(s) || val == __skip_value) continue;
            int spl = pl.create_edges(1);
            out_attr[spl] = attr[s];
            FOR(lv, 2) pl.vert(spl, lv) = m.vert(s, lv);
        }

        DropPolyLine(pl)._just_save_filename(__just_save_filename).add(out_attr, "out_attr")._just_save_filename(__just_save_filename)._show_vertices(false)._active_segment_attribute("out_attr").apply(name);
    }
};


template<PolyLineConceptBase C>
struct Drop<C, EdgeAttribute<vec3> > {
    Drop(C& pl, EdgeAttribute<vec3>& attr) : pl(pl), attr(attr), __forced_length(0), __forced_radius(-1) {}
    C& pl;
    EdgeAttribute<vec3>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(bool, use_disk, false);
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_length(double length, bool relative_to_bbox = false) {
        __forced_length = length;
        if (relative_to_bbox) {
            auto bbox = ToolBox(pl.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_length *= size;
        }
        return *this;
    };
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(pl.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    double __forced_length;
    void apply_arrow(std::string name = "") {
        apply(name, false);
    }
    void apply_disk(std::string name = "") {
        apply(name, true);
    }
protected:
    void apply(std::string name = "", bool use_disk = false) {
        PointSet middle;
        PointAttribute<vec3> middle_attr(middle);
        middle.create_points(pl.nedges());
        FOR(s, pl.nedges()) middle[s] = .5 * (pl.points[pl.vert(s, 0)] + pl.points[pl.vert(s, 1)]);
        FOR(s, pl.nedges()) middle_attr[s] = attr[s];
        Drop<PointSet, PointAttribute<vec3> >(middle, middle_attr)._just_save_filename(__just_save_filename)._force_length(__forced_length, false)._force_radius(__forced_radius, false)._resolution(__resolution).apply(name, use_disk);
    }
};




//     _____                 __                       
//    / ____|               / _|                      
//   | (___  _   _   _ __  | |_    __ _    ___    ___ 
//    \___ \| | | | | '__| |  _|  / _` |  / __|  / _ \
//    ____) | |_| | | |    | |   | (_| | | (__  |  __/
//   |_____/ \__,_| |_|    |_|    \__,_|  \___|  \___|
//                                                    
//                                                    

template <typename T> concept SurfaceConceptBase = std::is_base_of<Surface, T>::value;

inline double surface_ave_edge_size(Surface& m) {
    double sum = 0;
    FOR(f, m.nfacets()) FOR(lv, m.facet_size(f)) sum += (m.points[m.vert(f, lv)] - m.points[m.vert(f, (lv + 1) % m.facet_size(f))]).norm();
    return sum / double(m.ncorners());
}


template<SurfaceConceptBase C, class T>
struct Drop<C, FacetAttribute<T>> {
    Drop(C& m, FacetAttribute<T>& attr) : m(m), attr(attr) {}
    C& m;
    FacetAttribute<T>& attr;
    ARG(T, skip_value, T(-1));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return true && __skip_value != attr[v]; }));
    ARG(bool, show_edges, false);
    ARG(std::string, just_save_filename, "");
    void apply(std::string name = "") {
        Polygons other;
        FacetAttribute<T> other_attr(other);
        ToolBox(other.points).copy_from(m.points);
        FOR(f, m.nfacets()) if (__select(f)) {
            int nf = other.create_facets(1, m.facet_size(f));
            other_attr[nf] = attr[f];
            FOR(fv, m.facet_size(f)) other.vert(nf, fv) = m.vert(f, fv);
        }
        other.delete_isolated_vertices();
        if (__just_save_filename.compare("") == 0)
            DropSurface(other)._just_save_filename(__just_save_filename).add(other_attr, "attr")._show_vertices(false)._show_edges(__show_edges)._active_facet_attribute("attr").apply(name);
        else
            DropSurface(other)._just_save_filename(__just_save_filename).add(other_attr, "attr")._show_vertices(false)._show_edges(__show_edges)._active_facet_attribute("attr").apply(name);
    }
};

template<SurfaceConceptBase C, class T>
struct Drop<C, PointAttribute<T>> {
    Drop(C& m, PointAttribute<T>& attr) : m(m), attr(attr) {}
    C& m;
    PointAttribute<T>& attr;
    ARG(T, skip_value, T(-1));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(bool, show_edges, false);
    ARG(std::string, just_save_filename, "");
    void apply(std::string name = "") {
        Polygons other;
        PointAttribute<T> other_attr(other);
        ToolBox(other.points).copy_from(m.points);
        for (int v = 0; v < m.nverts(); v++) other_attr[v] = attr[v];
        FacetAttribute<bool> selected(m, true);
        for (auto h : m.iter_halfedges()) if (!__select(h.from()) || __skip_value == attr[h.from()]) selected[h.facet()] = false;

        FOR(f, m.nfacets()) if (selected[f]) {
            int nf = other.create_facets(1, m.facet_size(f));
            FOR(fv, m.facet_size(f)) other.vert(nf, fv) = m.vert(f, fv);
        }
        other.delete_isolated_vertices();
        DropSurface(other)._just_save_filename(__just_save_filename).add(other_attr, "attr")._show_vertices(false)._show_edges(__show_edges)._active_facet_attribute("attr").apply(name);
    }
};

template<SurfaceConceptBase C, class T>
struct Drop<C, CornerAttribute<T>> {
    Drop(C& m, CornerAttribute<T>& attr) : m(m), attr(attr), __forced_radius(-1) {}
    C& m;
    CornerAttribute<T>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(double, decal, .1);
    ARG(T, skip_value, T(-1));
    ARG(int, resolution, 10);
    ARG(bool, wireframe, false);
    ARG(std::string, just_save_filename, "");


    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;

    ARG(int, render_style, 0); // see code below 
    void apply(std::string name = "") {
        switch (__render_style) {
        case 0: apply_half_edge(name); break;
        case 1: apply_dual_edge(name); break;
        case 2: apply_corner(name); break;
        case 3: apply_interpolated(name); break;
        }
    }

    void apply_half_edge(std::string name = "") {
        um_assert(m.connected());
        PolyLine pl;
        EdgeAttribute<double> 			 out_attr(pl);


        for (auto h : m.iter_halfedges()) {
            vec3 G = Poly3(h.facet()).bary_verts();
            if (!__select(h) || __skip_value == attr[h]) continue;
            vec3 A = (1. - __decal) * m.points[h.from()] + __decal * G;
            vec3 B = (1. - __decal) * m.points[h.to()] + __decal * G;

            int e = pl.create_edges(1);
            out_attr[e] = attr[h];
            int v = pl.points.create_points(2);
            FOR(i, 2) pl.vert(e, i) = v + i;
            pl.points[v] = A;
            pl.points[v + 1] = B;
        }

        if (__wireframe)
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value).apply_wireframe(name);
        else
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._resolution(__resolution)._skip_value(__skip_value)._force_radius(__forced_radius).apply_arrow(name);
    }
    void apply_edge(std::string name = "") {
        um_assert(m.connected());
        PolyLine pl;
        EdgeAttribute<double> 			 out_attr(pl);


        for (auto h : m.iter_halfedges()) {
            auto opp = h.opposite();
            if (opp.active() && opp < h) continue;

            vec3 G = Poly3(h.facet()).bary_verts();
            if (!__select(h) || __skip_value == attr[h]) continue;
            vec3 A = m.points[h.from()];
            vec3 B = m.points[h.to()];

            int e = pl.create_edges(1);
            out_attr[e] = attr[h];
            int v = pl.points.create_points(2);
            FOR(i, 2) pl.vert(e, i) = v + i;
            pl.points[v] = A;
            pl.points[v + 1] = B;
        }

        if (__wireframe)
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value).apply_wireframe(name);
        else
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._resolution(__resolution)._skip_value(__skip_value)._force_radius(__forced_radius).apply_arrow(name);
    }

    void apply_dual_halfedge(std::string name = "") {
        um_assert(m.connected());
        PolyLine pl;
        EdgeAttribute<double> 			 out_attr(pl);
        for (auto h : m.iter_halfedges()) {
            vec3 G = Poly3(h.facet()).bary_verts();
            vec3 Gopp = .5 * (m.points[h.from()] + m.points[h.to()]);
            auto opp = h.opposite();
            if (opp.active()) Gopp = Poly3(opp.facet()).bary_verts();

            if (!__select(h) || __skip_value == attr[h]) continue;
            vec3 A = __decal * m.points[h.to()] + (1. - __decal) * G;
            vec3 B = __decal * m.points[h.to()] + (1. - __decal) * Gopp;
            int e = pl.create_edges(1);
            out_attr[e] = attr[h];
            int v = pl.points.create_points(2);
            FOR(i, 2) pl.vert(e, i) = v + i;
            pl.points[v] = A;
            pl.points[v + 1] = B;
        }
        if (__wireframe)
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value).apply_wireframe(name);
        else
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._resolution(__resolution)._skip_value(__skip_value)._force_radius(__forced_radius).apply_arrow(name);
        //Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value)._resolution(__resolution)._force_radius(__forced_radius).apply_arrow(name);
    }
    void apply_dual_edge(std::string name = "") {
        um_assert(m.connected());
        PolyLine pl;
        EdgeAttribute<double> 			 out_attr(pl);
        for (auto h : m.iter_halfedges()) {
            vec3 G = Poly3(h.facet()).bary_verts();
            vec3 Gopp = .5 * (m.points[h.from()] + m.points[h.to()]);
            auto opp = h.opposite();
            if (opp.active()) {
                Gopp = Poly3(opp.facet()).bary_verts();
                if (opp > h) continue;
            }
            if (!__select(h) || __skip_value == attr[h]) continue;
            vec3 A = G;
            vec3 B = Gopp;
            int e = pl.create_edges(1);
            out_attr[e] = attr[h];
            int v = pl.points.create_points(2);
            FOR(i, 2) pl.vert(e, i) = v + i;
            pl.points[v] = A;
            pl.points[v + 1] = B;
        }
        if (__wireframe)
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value).apply_wireframe(name);
        else
            Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._resolution(__resolution)._skip_value(__skip_value)._force_radius(__forced_radius).apply_arrow(name);
        //Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value)._resolution(__resolution)._force_radius(__forced_radius).apply_arrow(name);
    }
    void apply_corner(std::string name = "") {
        Triangles outm;
        FacetAttribute<T> out_attr(outm);
        um_assert(m.connected());
        for (auto h : m.iter_halfedges()) {
            T val = attr[h];
            if (!__select(h) || __skip_value == attr[h]) continue;
            vec3 pts[3] = {
                m.points[m.facets[h]],
                m.points[m.facets[h.next()]],
                m.points[m.facets[h.prev()]]
            };
            int offv = outm.points.create_points(3);
            outm.points[offv] = pts[0];
            outm.points[offv + 1] = 0.7 * pts[0] + .3 * pts[1];
            outm.points[offv + 2] = 0.7 * pts[0] + .3 * pts[2];
            int f = outm.create_facets(1);
            FOR(i, 3) outm.facets[3 * f + i] = offv + i;
            out_attr[f] = val;
        }
        DropSurface(outm)._just_save_filename(__just_save_filename).add(out_attr, "out_attr")._lighting(false)._show_border(false)._show_vertices(false)._active_facet_attribute("out_attr").apply(name);
    }
    void apply_interpolated(std::string name = "") {
        Polygons other;
        CornerAttribute<T> other_attr(other);
        ToolBox(other.points).copy_from(m.points);
        FacetAttribute<bool> selected(m, true);
        for (auto h : m.iter_halfedges()) if (!__select(h) || __skip_value == attr[h]) selected[h.facet()] = false;

        FOR(f, m.nfacets()) if (selected[f]) {
            int nf = other.create_facets(1, m.facet_size(f));

            FOR(fv, m.facet_size(f)) other_attr[other.corner(nf, fv)] = attr[m.corner(f, fv)];
            FOR(fv, m.facet_size(f)) other.vert(nf, fv) = m.vert(f, fv);
        }
        other.delete_isolated_vertices();
        DropSurface(other)._just_save_filename(__just_save_filename).add(other_attr, "attr")._show_vertices(false)._show_edges(false)._active_facet_attribute("attr").apply(name);
    }
};





template<SurfaceConceptBase C> struct Drop<C, CornerAttribute<vec3> > {
    Drop(C& m, CornerAttribute<vec3>& attr) : m(m), attr(attr), __forced_length(0), __forced_radius(-1) {}
    C& m;
    CornerAttribute<vec3>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(bool, use_disk, false);
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_length(double length, bool relative_to_bbox = false) {
        __forced_length = length;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_length *= size;
        }
        return *this;
    };
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    double __forced_length;
    void apply_arrow(std::string name = "") {
        apply(name, false);
    }
    void apply_disk(std::string name = "") {
        apply(name, true);
    }
protected:
    void apply(std::string name = "", bool use_disk = false) {
        PointSet pts;
        PointAttribute<vec3> pts_attr(pts);

        for (auto h : m.iter_halfedges())if (__select(h)) {
            int v = pts.create_points(1);
            pts[v] = .5 * (h.from().pos() + h.to().pos());
            pts_attr[v] = attr[h];
        }
        Drop<PointSet, PointAttribute<vec3> >(pts, pts_attr)._just_save_filename(__just_save_filename)._force_length(__forced_length, false)._force_radius(__forced_radius, false)._resolution(__resolution).apply(name, use_disk);
    }
};


template<SurfaceConceptBase C> struct Drop<C, FacetAttribute<vec3> > {
    Drop(C& m, FacetAttribute<vec3>& attr) : m(m), attr(attr), __forced_length(0), __forced_radius(-1) {}
    C& m;
    FacetAttribute<vec3>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(bool, use_disk, false);
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_length(double length, bool relative_to_bbox = false) {
        __forced_length = length;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_length *= size;
        }
        return *this;
    };
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    double __forced_length;
    void apply_arrow(std::string name = "") {
        apply(name, false);
    }
    void apply_disk(std::string name = "") {
        apply(name, true);
    }
    void apply_line(std::string name = "") {
        apply(name, false, true);
    }
protected:
    void apply(std::string name = "", bool use_disk = false, bool use_line = true) {
        PointSet pts;
        PointAttribute<vec3> pts_attr(pts);

        for (auto f : m.iter_facets())if (__select(f)) {
            int v = pts.create_points(1);
            pts[v] = Poly3(f).bary_verts();
            pts_attr[v] = attr[f];
        }
        if (use_line)
            Drop<PointSet, PointAttribute<vec3> >(pts, pts_attr)._just_save_filename(__just_save_filename)._force_length(__forced_length, false)._force_radius(__forced_radius, false)._resolution(__resolution).apply(name, use_disk, use_line);
    }

};







template<SurfaceConceptBase C> struct Drop<C, CornerAttribute<vec2> > {
    Drop(C& m, CornerAttribute<vec2>& attr) : m(m), attr(attr), __forced_radius(-1) {}
    C& m;
    CornerAttribute<vec2>& attr;
    ARG(double, radius_wrt_bbox, 0);
    ARG(std::string, just_save_filename, "");
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    void apply_texture(std::string name = "") {
        DropSurface(m)._just_save_filename(__just_save_filename).add(attr, "attr")._show_vertices(false)._active_uv_attribute("attr").apply(name);
    }
    void apply_neg_detJ(std::string name = "") {
        FacetAttribute<bool> has_neg_J(m);
        FOR(f, m.nfacets()) has_neg_J[f] = (mat2x2({ attr[m.corner(f, 1)] - attr[m.corner(f, 0)], attr[m.corner(f, 2)] - attr[m.corner(f, 0)] }).det() <= 0);
        int nb_neg = 0;
        FOR(f, m.nfacets()) if (has_neg_J[f]) nb_neg++;


        Trace::log_value("det_J_inf_0_" + name, double(nb_neg) / double(m.nfacets()));

        if (nb_neg == 0) return;
        Drop<C, FacetAttribute<bool> >(m, has_neg_J)._just_save_filename(__just_save_filename)._skip_value(false).apply(name);
    }
    void apply_map(std::string name = "") {
        Triangles out;
        out.points.create_points(3 * m.nfacets());
        out.create_facets(m.nfacets());
        FOR(v, 3 * m.nfacets()) {
            out.vert(v / 3, v % 3) = v;
            out.points[v] = vec3(attr[v][0], attr[v][1], 0);
        }
        DropSurface(out)._just_save_filename(__just_save_filename).apply(name);
    }
    void apply_iso_int(std::string name = "") {

        um_assert(m.connected());
        PolyLine poly;
        EdgeAttribute<int> iso_coord(poly);
        for (auto d : { 0,1 }) {
            for (auto f : m.iter_facets()) {
                //  plop(f);
                double minv = 1e20;
                double maxv = -1e20;
                for (auto cir : f.iter_halfedges()) {
                    minv = std::min(minv, attr[cir][d]);
                    maxv = std::max(maxv, attr[cir][d]);
                }
                for (double iso = std::floor(minv); iso < std::ceil(maxv) + 1; iso += 1.) {
                    std::vector<vec3> pts;

                    for (auto cir : f.iter_halfedges()) {
                        if (std::abs(attr[cir.next()][d] - attr[cir][d]) < 1e-3)continue;
                        double c = (iso - attr[cir][d]) / (attr[cir.next()][d] - attr[cir][d]);
                        if (c < -.0001) continue;
                        if (c > 1.0001) continue;
                        pts.push_back((1. - c) * m.points[cir.from()] + c * m.points[cir.to()]);
                    }
                    if (pts.size() == 3) {
                        if ((pts[0] - pts[1]).norm2() < 1e-5) std::swap(pts[0], pts[2]);
                        pts.pop_back();
                    }
                    if (pts.size() == 2) {
                        int offs = poly.create_edges(1);
                        iso_coord[offs] = d;
                        int offv = poly.points.create_points(2);
                        FOR(e, 2) poly.vert(offs, e) = offv + e;
                        FOR(e, 2) poly.points[offv + e] = pts[e];
                    }
                }
            }
        }
        Drop<PolyLine, EdgeAttribute<int>>(poly, iso_coord)._just_save_filename(__just_save_filename)._force_radius(__forced_radius, false).apply_capsule(name);
    }
};





template<SurfaceConceptBase C> struct Drop<C, PointAttribute<vec2> > {
    Drop(C& m, PointAttribute<vec2>& attr) : m(m), attr(attr), __forced_radius(-1) {}
    C& m;
    PointAttribute<vec2>& attr;
    ARG(double, radius_wrt_bbox, 0);
    ARG(std::string, just_save_filename, "");
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    void apply_texture(std::string name = "") {
        DropSurface(m)._just_save_filename(__just_save_filename).add(attr, "attr")._active_point_attribute("attr").apply(name);
    }
};





template<SurfaceConceptBase C> struct Drop<C, PointAttribute<vec3> > {
    Drop(C& m, PointAttribute<vec3>& attr) : m(m), attr(attr){}
    C& m;
    PointAttribute<vec3>& attr;
    ARG(std::string, just_save_filename, "");
    void apply(std::string name = "") {
        DropSurface(m)._active_point_attribute("attr[0]").add(attr, "attr")._just_save_filename(__just_save_filename).apply(name);
    }
};
//   __      __     _                            
//   \ \    / /    | |                           
//    \ \  / /__   | |  _   _   _ __ ___     ___ 
//     \ \/ / _ \  | | | | | | | '_ ` _ \   / _ \
//      \  / (_) | | | | |_| | | | | | | | |  __/
//       \/ \___/  |_|  \__,_| |_| |_| |_|  \___|
//                                               
//                                              




template <typename T> concept VolumeConceptBase = std::is_base_of<Volume, T>::value;

template<VolumeConceptBase C, class  T>
struct Drop<C, CellAttribute<T> > {
    Drop(C& m, CellAttribute<T>& attr) : m(m), attr(attr) {}
    C& m;
    CellAttribute<T>& attr;
    ARG(T, skip_value, T(-1));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(std::string, just_save_filename, "");
    void apply(std::string name = "") {
        C cpy;
        ToolBox(cpy).copy_from(m);

        CellAttribute<T> cpy_attr(cpy);
        FOR(c, m.ncells())	cpy_attr[c] = attr[c];
        std::vector<bool>	to_kill(m.ncells());
        FOR(c, m.ncells())	to_kill[c] = (!__select(c) || attr[c] == __skip_value);
        cpy.delete_cells(to_kill);
        cpy.delete_isolated_vertices();
        if (cpy.ncells() == 0) std::cerr << "Output " << name << " has no cells, so it is skipped\n";
        else DropVolume(cpy)._just_save_filename(__just_save_filename).add(cpy_attr, "attr")._active_cell_attribute("attr").apply(name);
    }
};


void show_iso(Tetrahedra& m, CellCornerAttribute<vec3>& U, std::string name = "iso", int nb_subdiv = 1);
void show_neg_det_J(Tetrahedra& m, CellCornerAttribute<vec3>& U, std::string name);
void show_U_deformation(Tetrahedra& tet, CellCornerAttribute<vec3>& U);

template<>
struct Drop<Tetrahedra, CellCornerAttribute<vec3> > {
    Drop(Tetrahedra& m, CellCornerAttribute<vec3>& U) : m(m), U(U) {}
    Tetrahedra& m;
    CellCornerAttribute<vec3>& U;
    ARG(vec3, skip_value, vec3(0, 0, 0));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(std::string, just_save_filename, "");
    Drop<Tetrahedra, CellCornerAttribute<vec3> >& apply(std::string name = "") {
        DropVolume(m).add(U, "U")._just_save_filename(__just_save_filename)._active_cell_corner_attribute("U[0]").apply(name);
        return *this;
    }
    Drop<Tetrahedra, CellCornerAttribute<vec3> >& apply_neg_det(std::string name = "negdet") {
        show_neg_det_J(m, U, name);
        return *this;
    }
    Drop<Tetrahedra, CellCornerAttribute<vec3> >& apply_detJ(std::string name = "detJ") {
        CellAttribute<double> detJ(m);
        for (auto c : m.iter_cells()) {
            mat3x3 J;
            Tetrahedron tet = Tetrahedron(c);
            mat<3, 4> grd = tet.grad_operator();
            FOR(d, 3) {
                vec4 scal;
                FOR(lv, 4) scal[lv] = U[c.corner(lv)][d];
                J[d] = grd * scal;
            }
            detJ[c] = J.det();
        }
        Drop<Tetrahedra, CellAttribute<double> >(m, detJ)._just_save_filename(__just_save_filename)._skip_value(1e20).apply(name);

        return *this;
    }
    Drop<Tetrahedra, CellCornerAttribute<vec3> >& apply_iso(std::string name = "", int frequency = 1) {
        show_iso(m, U, name, frequency);
        return *this;
    }



    Drop<Tetrahedra, CellCornerAttribute<vec3> >& apply_deformation(std::string name = "") {
        show_U_deformation(m, U);
        return *this;
    }



    Drop<Tetrahedra, CellCornerAttribute<vec3> >& apply_map(std::string name = "") {

        Tetrahedra tet;
        tet.points.create_points(4 * m.ncells());
        tet.create_cells(m.ncells());
        for (auto c : m.iter_cells()) FOR(lv, 4) {
            tet.vert(c, lv) = 4 * c + lv;
            tet.points[4 * c + lv] = U[c.corner(lv)];
        }
        DropVolume(tet)._just_save_filename(__just_save_filename).apply(name);
        return *this;
    }
};




template<VolumeConceptBase C>
struct Drop<C, CellAttribute<vec3>> {
    Drop(C& m, CellAttribute<vec3>& attr) : m(m), attr(attr), __forced_length(0), __forced_radius(-1) {}
    C& m;
    CellAttribute<vec3>& attr;
    ARG(std::function< bool(int) >, select, std::function([&](int c) {return attr[c].norm2() != 0; }));
    ARG(bool, use_disk, false);
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_length(double length, bool relative_to_bbox = false) {
        __forced_length = length;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_length *= size;
        }
        return *this;
    };
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    double __forced_length;
    void apply_arrow(std::string name = "") {
        apply(name, false);
    }
    void apply_disk(std::string name = "") {
        apply(name, true);
    }
protected:
    void apply(std::string name = "", bool use_disk = false) {
        PointSet pts;
        PointAttribute<vec3> pts_attr(pts);

        for (Volume::Cell c : m.iter_cells()) if (__select(c)) {
            int v = pts.create_points(1);
            pts[v] = vec3(0, 0, 0);
            FOR(lv, c.nverts())pts[v] += (1. / double(c.nverts())) * c.vertex(lv).pos();
            pts_attr[v] = attr[c];
        }
        Drop<PointSet, PointAttribute<vec3> >(pts, pts_attr)._just_save_filename(__just_save_filename)._force_length(__forced_length, false)._force_radius(__forced_radius, false)._resolution(__resolution).apply(name, use_disk);
    }

};







template<VolumeConceptBase C, class T>
struct Drop<C, CellCornerAttribute<T>> {
    Drop(C& m, CellCornerAttribute<T>& attr) : m(m), attr(attr) {}
    C& m;
    CellCornerAttribute<T>& attr;
    ARG(std::string, just_save_filename, "");
    void apply(std::string name = "") {
        DropVolume(m).add(attr, "attr")._just_save_filename(__just_save_filename)._active_cell_corner_attribute("attr").apply(name);
    }
};

template<VolumeConceptBase C, class T>
struct Drop<C, PointAttribute<T>> {
    Drop(C& m, PointAttribute<T>& attr) : m(m), attr(attr) {}
    C& m;
    PointAttribute<T>& attr;
    ARG(std::string, just_save_filename, "");
    void apply(std::string name = "") {
        DropVolume(m).add(attr, "attr")._just_save_filename(__just_save_filename)._active_point_attribute("attr").apply(name);
    }
    void apply_iso(std::string name = "",int nb_iso=10) {
        Triangles iso;
        FacetAttribute<double> surf_attr(iso);
        CellCornerAttribute<double> scal(m);

        double range[2] = { 1e20,-1e20 };
        for (auto v : m.iter_vertices()) {
            range[0] = std::min(range[0], attr[v]);
            range[1] = std::max(range[1], attr[v]);
        }

        for (auto c : m.iter_corners()) scal[c] = -.001 + 1.02 * nb_iso * (attr[c.vertex()] - range[0]) / (range[1] - range[0]);
        FOR(i, nb_iso + 1) extract_iso(m, iso, surf_attr, scal, i);
        Drop<Surface, FacetAttribute<double>>(iso, surf_attr).apply(name);
    }

};


template<VolumeConceptBase C, class T>
struct Drop<C, CellFacetAttribute<T>> {
    Drop(C& m, CellFacetAttribute<T>& attr) : m(m), attr(attr), __forced_radius(-1) {}
    C& m;
    CellFacetAttribute<T>& attr;
    ARG(double, shrink, 0);
    ARG(T, skip_value, T(-1));
    ARG(std::function< bool(int) >, select, std::function([&](int v) {return v!=-1; }));
    ARG(int, resolution, 10);
    ARG(std::string, just_save_filename, "");
    auto& _force_radius(double radius, bool relative_to_bbox = false) {
        __forced_radius = radius;
        if (relative_to_bbox) {
            auto bbox = ToolBox(m.points).bbox();
            double size = (bbox.max - bbox.min).norm();
            __forced_radius *= size;
        }
        return *this;
    };
    double __forced_radius;
    void apply(std::string name = "") {
        Polygons surf;
        FacetAttribute<T> surf_attr(surf);
        //OppositeFacet adjacent(m);
        for (auto c : m.iter_cells()) FOR(cf, m.nfacets_per_cell()) {
            if (!__select(m.facet(c, cf)) || __skip_value == attr[m.facet(c, cf)]) continue;
            int off_v = surf.points.create_points(m.facet_size(cf));
            vec3 G = vec3(0, 0, 0);
            FOR(lv, c.nverts()) G += (1. / double(c.nverts())) * c.vertex(lv).pos();
            FOR(cfv, m.facet_size(cf)) surf.points[off_v + cfv] = __shrink * G + (1. - __shrink) * m.points[m.facet_vert(c, cf, cfv)];
            int off_f = surf.create_facets(1, m.facet_size(cf));
            FOR(cfv, m.facet_size(cf)) surf.vert(off_f, cfv) = off_v + cfv;
            surf_attr[off_f] = attr[m.facet(c, cf)];
        }

        Drop<Polygons, FacetAttribute<T>>(surf, surf_attr)._just_save_filename(__just_save_filename)._skip_value(__skip_value).apply(name);
    }
    void apply_dual(std::string name = "") {

        PolyLine pl;
        EdgeAttribute<double> 			 out_attr(pl);
        FOR(c, m.ncells()) FOR(cf, m.nfacets_per_cell()) {
            vec3 G_c;
            vec3 G_f;
            if (std::is_same<C, Hexahedra>::value) {
                G_c = Hexahedron(Volume::Cell(m, c)).bary_verts();
                G_f = Quad3(Volume::Cell(m, c).facet(cf)).bary_verts();
            }
            else if (std::is_same<C, Tetrahedra>::value) {
                G_c = Tetrahedron(Volume::Cell(m, c)).bary_verts();
                G_f = Triangle3(Volume::Cell(m, c).facet(cf)).bary_verts();
            }
            else um_assert(!"unknowned type of volume mesh");

            if (__skip_value != attr[m.facet(c, cf)]) {
                int e = pl.create_edges(1);
                out_attr[e] = attr[m.facet(c, cf)];
                int v = pl.points.create_points(2);
                FOR(i, 2) pl.vert(e, i) = v + i;
                pl.points[v] = G_c;//m.util.bary_verts(c);
                pl.points[v + 1] = G_f;//m.util.bary_facet(c, cf);
            }
        }
        Drop<PolyLine, EdgeAttribute<double>>(pl, out_attr)._just_save_filename(__just_save_filename)._resolution(__resolution)._force_radius(__forced_radius).apply_capsule(name);
    }
};


template<VolumeConceptBase C>
struct Drop<C, CellAttribute<mat3x3>> {
    Drop(C& m, CellAttribute<mat3x3>& J) : m(m), J(J) {}
    C& m;
    CellAttribute<mat3x3>& J;
    ARG(std::function< bool(int) >, select, std::function([&](int c) {return J[c].norm() != 0; }));
    ARG(double,scale,1);
    ARG(bool,with_flag,false);
    ARG(std::string, just_save_filename, "");
    void apply_cube(std::string name = "") {
          double ave = ToolBox(m).ave_edge_size();
        Hexahedra outm;
        outm.create_cells(m.ncells());
        outm.points.create_points(8 * m.ncells());
        CellFacetAttribute<int> flag(outm, 0);
        FOR(f, outm.nfacets())flag[f] = f % 6;
    for (auto c : m.iter_cells()) {
        FOR(v, 8) outm.vert(c, v) = 8 * c + v;
        FOR(di, 2)FOR(dj, 2)FOR(dk, 2) outm.points[8 * c + di + 2 * dj + 4 * dk] 
            = Tetrahedron(c).bary_verts() + __scale * .15 * ave
            * J[c].transpose() /*J is a rot => transpose == fast_invert*/ * vec3(2 * di - 1, 2 * dj - 1, 2 * dk - 1);
    }
    if (__with_flag)	Drop<Hexahedra,CellFacetAttribute<int> >(outm, flag).apply(name + "_cubes");
    else			    DropVolume(outm).apply(name + "_cubes");

    }

};





#endif