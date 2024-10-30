#ifndef __META_MESH__H__
#define __META_MESH__H__

#include <filesystem>
#include <framework/trace.h>
#include <fullhex/cube_cover.h>
#include <fullhex/fiber.h>
#include <fullhex/frame_field_3d.h>
#include <fullhex/gp_basic.h>
#include <fullhex/hextract.h>
#include <volume/frame3D.h>
#include <volume/hex_stat.h>
#include <create_shape.h>

namespace Meta {
    template<class C>
    struct Mesh2d {
        Mesh2d(C& m, std::string filename = ""):m(m) {
            if (filename.empty() ) return;
            if (!std::filesystem::exists(filename)) return;
            attr = read_by_extension(filename, m);
            m.connect();
            Trace::log_value("#facets", m.nfacets());
            Trace::log_value("#vert", m.nverts());
        }
        bool active() { return m.nverts() > 0; }

        template<class T> Mesh2d& add(PointAttribute<T>& c, std::string str) { attr.points.push_back({ str,c.ptr }); return *this; }
        template<class T> Mesh2d& add(CornerAttribute<T>& c, std::string str) { attr.corners.push_back({ str,c.ptr }); return *this; }
        template<class T> Mesh2d& add(FacetAttribute<T>& c, std::string str) { attr.facets.push_back({ str,c.ptr }); return *this; }

        template<class S> Mesh2d& get_copy(S& c, std::string str) {
            S tmp(str, attr, m); std::copy(tmp.ptr->data.begin(), tmp.ptr->data.end(), c.ptr->data.begin());
            return *this;
        }
        template<class S> Mesh2d& get(S& c, std::string str) {
            S tmp(str, attr, m);
            c.ptr = tmp.ptr;
            return *this;
        }

        void save(std::string filename) { write_by_extension(filename, m, attr); }

        SurfaceAttributes attr = { {},{},{} };
        C& m;
    };


    template<class C>
    struct Mesh3d {
        Mesh3d(C& m,std::string filename = "") : m(m) {
            if (!filename.empty()) {
                attr = read_by_extension(filename, m);
                FOR(i, attr.cell_corners.size()) {
                    std::string full_attr_name = attr.cell_corners[i].first;
                    if (full_attr_name.length() < 8) continue;
                    std::string attr_name = full_attr_name.substr(0, full_attr_name.length() - 8);
                    if (full_attr_name == attr_name + "_as_vec3") {

                        std::cerr   <<"Import mat3x3 cell attribute "<< attr_name
                                    <<"  From vec3 cellcorner attribute "<< full_attr_name<<std::endl;
                        CellAttribute<mat3x3> J(attr_name, attr, m);
                        CellCornerAttribute<vec3> tmp(full_attr_name, attr, m);
                        for (auto c : m.iter_cells()) J[c] = uvw_to_jacobian(c, tmp);
                        std::swap(attr.cell_corners[i], attr.cell_corners.back());
                        attr.cell_corners.pop_back();
                    }

                }
                m.connect();
            }
        }
        bool active() { return m.nverts() > 0; }

        template<class T> Mesh3d& add(PointAttribute<T>& c, std::string str) { 
            for (auto a:attr.points) if (a.first==str) {
                std::cerr<<str<<" is already a registered attribute\n";
                return *this;
            }
            attr.points.push_back({ str,c.ptr }); 
            return *this; 
        }
        template<class T> Mesh3d& add(CellAttribute<T>& c, std::string str) {
            for (auto a:attr.cells) if (a.first==str) {
                std::cerr<<str<<" is already a registered attribute\n";
                return *this;
            }
            attr.cells.push_back({ str,c.ptr });
            return *this;
        }
        template<class T> Mesh3d& add(CellFacetAttribute<T>& c, std::string str) { 
            for (auto a:attr.cell_facets) if (a.first==str) {
                std::cerr<<str<<" is already a registered attribute\n";
                return *this;
            }
            attr.cell_facets.push_back({ str,c.ptr }); 
            return *this; 
        }
        template<class T> Mesh3d& add(CellCornerAttribute<T>& c, std::string str) { 
            for (auto a:attr.cell_corners) if (a.first==str) {
                std::cerr<<str<<" is already a registered attribute\n";
                return *this;
            }
            attr.cell_corners.push_back({ str,c.ptr }); 
            return *this; 
        }

        template<class S> Mesh3d& get_copy(S& c, std::string str) {
            S tmp(str, attr, m); std::copy(tmp.ptr->data.begin(), tmp.ptr->data.end(), c.ptr->data.begin());
            return *this;
        }

        inline void convert_mat3x3_to_vec3() {
            int i = 0;
            while (i < attr.cells.size()) {
                auto ptr = dynamic_cast<std::shared_ptr<AttributeContainer<mat3x3>>::element_type*> (attr.cells[i].second.get());
                if (ptr == NULL) { i++; continue; }
                std::string mat_name = attr.cells[i].first;

                std::cerr << " ptr " << ptr << " attr names " << attr.cells[i].first << std::endl;

                CellAttribute<mat3x3> J(mat_name, attr, m);
                CellCornerAttribute<vec3> tmp(mat_name + "_as_vec3",attr,m);
                for (auto c : m.iter_cells()) {
                    mat<4, 3> val = jacobian_to_uvw(c, J);
                    FOR(lv, 4) tmp[c.corner(lv)] = val[lv];
                }
                std::swap(attr.cells[i], attr.cells.back());
                attr.cells.pop_back();
                i = 0;
            }
        }

        inline void save(std::string filename) {
            convert_mat3x3_to_vec3();
             for (auto it : attr.cells) plop(it.first);
            write_by_extension(filename, m, attr);
        }


        inline void drop(std::string name) {
            convert_mat3x3_to_vec3();

           if (!Trace::drop_mesh_is_active) return ;
		    std::string filename = Trace::outputdir + "/" + name + "_" + std::to_string(Trace::num_drop++) + ".geogram";
		    write_by_extension(filename, m, attr);
		std::ofstream myfile;
		myfile.open(Trace::outputdir + "/view.lua", std::ios::out | std::ios::app);
		myfile << "scene_graph.load_object(\"" << filename << "\")\n";
		myfile << "scene_graph.current().shader.vertices_style = 'true; 0 0 0 1; 1'\n";
		myfile << "scene_graph.current().shader.shrink = 1\n";
		myfile.close();
        }
        VolumeAttributes attr = {};
        C &m;
    };
};




inline void compute_fid(Meta::Mesh2d<Triangles>& mesh) {
    Trace::step("compute_fid");
    FacetAttribute<int> fid(mesh.m);
    std::iota(fid.ptr->data.begin(), fid.ptr->data.end(), 0);
    mesh.add(fid, "fid");
}


inline void sum_fid_in_one_ring(Meta::Mesh2d<Triangles>& mesh) {
    Trace::step("sum_fid_in_one_ring");
    // unpack input mesh2d
    Triangles& m = mesh.m;				m.connect();
    FacetAttribute<int> fid(m);			mesh.get(fid, "fid");

    // add output
    PointAttribute<double> ave(m, 0);	mesh.add(ave, "ave");

    // do some stuff
    for (auto v : m.iter_vertices()) for (auto cir : v.iter_halfedges()) ave[v] += fid[cir.facet()];
    Drop(m, ave).apply();
}




inline void test_surface() {
    std::string filename = "C:\\NICO\\data\\catorus.geogram";
    std::string core = "C:\\NICO\\tmp\\last.geogram";


    Triangles m;
    if (0)
    {
        Meta::Mesh2d tri(m, core);
        //compute_fid(tri);
        sum_fid_in_one_ring(tri);
        return;
    }


    Meta::Mesh2d tri(m, filename);

    //tri.save(core);return;
    compute_fid(tri);
    FOR(i, 5) {
        //tri.save(core);return;
        sum_fid_in_one_ring(tri);
    }
}


inline void compute_cid(Meta::Mesh3d<Tetrahedra>& mesh) {
    Trace::step("compute_fid");
    CellAttribute<int> cid("cid",mesh.attr,mesh.m);
    std::iota(cid.ptr->data.begin(), cid.ptr->data.end(), 0);
    Drop(mesh.m, cid).apply();
    //mesh.add(cid, "cid");
}

inline void compute_U(Meta::Mesh3d<Tetrahedra>& mesh) {
    Tetrahedra& m = mesh.m;
    Trace::step("compute_fid");
    CellCornerAttribute<vec3> U("U",mesh.attr,mesh.m);
    for (auto c : m.iter_corners())
        U[c] = c.vertex().pos();
    Drop(m, U).apply_iso();

    CellAttribute<mat3x3> J("J",mesh.attr,mesh.m);
    for (auto c : m.iter_cells()) J[c] = mat3x3::identity();

    //mesh.add(J, "J");
}

inline void show_U(Meta::Mesh3d<Tetrahedra>& mesh) {
    Tetrahedra& m = mesh.m;
    Trace::step("show U");
    CellCornerAttribute<vec3> U("U",mesh.attr,mesh.m);
    for (auto c : m.iter_corners())
        U[c] = 2.*c.vertex().pos();
    Drop(m, U).apply_iso();
//    mesh.add(U, "U");

    CellAttribute<mat3x3> J("J", mesh.attr, mesh.m);
    Drop(mesh.m, J).apply_cube("J");
}



inline void test_volume() {
    std::string filename = "C:\\NICO\\data\\mesh3d-tet\\pente.mesh";
    std::string core = "C:\\NICO\\tmp\\last.geogram";

    if (0)
    {
        Tetrahedra tet;
        Meta::Mesh3d mesh(tet,core);
        compute_cid(mesh);
        return;
    }


    Tetrahedra tet;
    Meta::Mesh3d mesh(tet,filename);
    compute_cid(mesh);
    compute_U(mesh);
    show_U(mesh);
    mesh.save(core);
    Tetrahedra tet_load;
    Meta::Mesh3d mesh_load(tet_load,core);
    CellCornerAttribute<vec3> U("U", mesh_load.attr, tet_load);
    Drop(mesh_load.m, U).apply("U");
    show_U(mesh_load);
    //tet_load.save("droptest.geogram");
    mesh_load.drop("droptest");
}
#endif