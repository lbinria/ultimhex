#include <drop_attribute.h>
#include <geom3.h>
#include <fullhex/gp_basic.h>



DropPolyLineGeometry::DropPolyLineGeometry(PolyLine& m) : m(m), __forced_radius(-1) {}
DropPolyLineGeometry& DropPolyLineGeometry::_force_radius(double radius, bool relative_to_bbox ) {
		__forced_radius = radius;
		if (relative_to_bbox) {
			auto bbox = ToolBox(m.points).bbox();
			double size = (bbox.max - bbox.min).norm();
			__forced_radius *= size;
		}
		return *this;
	};

	void DropPolyLineGeometry::auto_set_radius() {
		if (__forced_radius != -1) return;
		vec2 ave_edge_length(0, 0);
		FOR(s, m.nedges())
			if ((m.points[m.vert(s, 0)] - m.points[m.vert(s, 1)]).norm2() > 1e-10)
				ave_edge_length += vec2((m.points[m.vert(s, 0)] - m.points[m.vert(s, 1)]).norm(), 1);
		__forced_radius = .01 * ave_edge_length[0] / ave_edge_length[1];
	}

	void DropPolyLineGeometry::apply_arrow(std::string name ) {
		Polygons outm;
		auto_set_radius();
		FOR(s, m.nedges())
			Glyphs::Arrow(outm).resolution(__resolution).radius(__forced_radius).apply(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)]);
		DropSurface(outm)._just_save_filename(__just_save_filename)._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false).apply(name);
	}

	void DropPolyLineGeometry::apply_capsule(std::string name) {
		Polygons outm;
		auto_set_radius();
		FOR(s, m.nedges())
			Glyphs::Capsule(outm).resolution(__resolution).radius(__forced_radius).apply(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)]);
		DropSurface(outm)._just_save_filename(__just_save_filename)._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false).apply(name);
	}
	void DropPolyLineGeometry::apply_disk(std::string name ) {
		Polygons outm;
		if (__forced_radius == -1) {
			auto_set_radius();
			__forced_radius *= 4.;
		}
		FOR(s, m.nedges())
			Glyphs::Disk(outm).resolution(__resolution).radius(__forced_radius).apply(m.points[m.vert(s, 0)], m.points[m.vert(s, 1)] - m.points[m.vert(s, 0)]);
		DropSurface(outm)._just_save_filename(__just_save_filename)._lighting(false)._show_border(false)._show_edges(false)._show_vertices(false).apply(name);
	}



void extract_iso(Tetrahedra& m, Triangles& tr, FacetAttribute<double>& attr, CellCornerAttribute<double>& scal, double iso) {

	//auto heh = m.heh;
	FOR(c, m.ncells()) //FOR(lf, 4) FOR(lv, 3)
	{
		vec3 P[4];
		double val[4];
		FOR(lv, 4) P[lv] =  m.points[m.vert(c,lv)];
		FOR(lv, 4) val[lv] =  scal[m.corner(c,lv)];
		FOR(lv, 4) val[lv] -= iso;

		auto swap_vert = [&](int i, int j){
			std::swap(val[i],val[j]);
			std::swap(P[i],P[j]);
		};

		{	
			bool allsup = true; 
			FOR(i, 4) allsup = (allsup && (val[i] >= 0)); 
			if (allsup) FOR(i, 4) val[i] -= 1e-10;
			else FOR(i, 4) val[i] += 1e-10;
		}
		// bubble sort val
		FOR(iter,3) FOR(lv,3)if (val[lv]>val[lv+1]) swap_vert(lv,lv+1);
		
		FOR(lv,3) um_assert(val[lv]<=val[lv+1]);

		if (val[0]>0 || val[3]<0) continue;

		if (val[1]*val[2]>0){
			if (val[2]<0) swap_vert(0,3);
			int offv = tr.points.create_points(3);
			FOR(i, 3) {
				int lv = i+1;
				double lambda = val[lv] / (val[lv] - val[0]);
				tr.points[offv + i] = lambda * P[0] + (1. - lambda) * P[lv];
			}
			int f = tr.create_facets(1);
			attr[f] = iso;
			FOR(i, 3) tr.vert(f, i) = offv + i;
		}

		if (val[2]>0 && val[1]<0) {
			int offv = tr.points.create_points(4);
			int edges[4][2] = { {0,3},{1,3},{1,2},{0,2} };
			FOR(e, 4) {
				double diff = val[edges[e][1]] - val[edges[e][0]];
				double lambda = val[edges[e][1]]  / diff;
				tr.points[offv + e] = lambda * P[edges[e][0]] + (1. - lambda) * P[edges[e][1]];
			}
			int f = tr.create_facets(2);
			FOR(i, 2) attr[f + i] = iso;
			tr.vert(f, 0) = offv;	tr.vert(f, 1) = offv + 1;		tr.vert(f, 2) = offv + 2;
			tr.vert(f + 1, 0) = offv; tr.vert(f + 1, 1) = offv + 2;	tr.vert(f + 1, 2) = offv + 3;
		}

		continue;



		//if (h[0] > h[1] && h[0] > h[2])
		{ // avoid duplicated triangles... duplicated vertices are fines;) 
			if (val[3] < 0 && val[0] > 0 && val[1] > 0 && val[2] > 0) {
				int offv = tr.points.create_points(3);
				FOR(i, 3) {
					double lambda = val[i] / (val[i] - val[3]);
					tr.points[offv + i] = lambda * P[3] + (1. - lambda) * P[i];
				}
				int f = tr.create_facets(1);
				attr[f] = iso;
				FOR(i, 3) tr.vert(f, i) = offv + i;
			}
			if (val[3] > 0 && val[0] < 0 && val[1] < 0 && val[2] < 0) {
				int offv = tr.points.create_points(3);
				FOR(i, 3) {
					double lambda = val[i] / (val[i] - val[3]);
					tr.points[offv + i] = lambda * P[3] + (1. - lambda) * P[i];
				}
				int f = tr.create_facets(1);
				attr[f] = iso;
				FOR(i, 3) tr.vert(f, 2 - i) = offv + i;
			}
		}

		if (val[3] > 0 && val[0] > 0 && val[1] < 0 && val[2] < 0) {
			int offv = tr.points.create_points(4);
			int edges[4][2] = { {0,1},{0,2},{3,2},{3,1} };
			FOR(e, 4) {
				double lambda = val[edges[e][1]]  / (val[edges[e][1]] - val[edges[e][0]]);
				tr.points[offv + e] = lambda * P[edges[e][0]] + (1. - lambda) * P[edges[e][1]];
			}
			int f = tr.create_facets(2);
			FOR(i, 2) attr[f + i] = iso;
			tr.vert(f, 0) = offv;	tr.vert(f, 1) = offv + 1;		tr.vert(f, 2) = offv + 2;
			tr.vert(f + 1, 0) = offv; tr.vert(f + 1, 1) = offv + 2;	tr.vert(f + 1, 2) = offv + 3;
		}
	}
}
template <class T> vec2 scalar_attribute_range(GenericAttribute<T>& attr) {
	vec2 res(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
	//vec2 res(1e20, -1e20);
	FOR(i, attr.ptr->data.size()) res[0] = std::min<double>(res[0], attr[i]);
	FOR(i, attr.ptr->data.size()) res[1] = std::max<double>(res[1], attr[i]);
	return res;
}

void show_iso(Tetrahedra& m,CellCornerAttribute<vec3>& U,std::string name,int nb_subdiv ){
	if (!Trace::drop_mesh_is_active) return;
		FOR(dim, 3) {
			CellCornerAttribute<double> scal(m);
			FOR(c, 4 * m.ncells()) scal[c] = U[c][dim]*double(nb_subdiv);
			vec2 range = scalar_attribute_range(scal);
			Triangles tr;
			FacetAttribute<double> attr(tr);
			for (int i = int(std::floor(range[0])); i < int(std::ceil(range[1])) + 1; i++)
				extract_iso(m, tr, attr, scal, i);
			Drop(tr, attr)._skip_value(1e20).apply(name);
		}
}

void show_neg_det_J(Tetrahedra& m, CellCornerAttribute<vec3>& U,std::string name) {
	if (!Trace::drop_mesh_is_active) return;
	CellAttribute<double> detJ(m);
	bool need_render = false;
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
		if (detJ[c] > 1e-5) detJ[c] = 1e20; 
		else need_render = true; 
	}
	if (need_render) Drop(m, detJ)._skip_value(1e20).apply(name);
}


void show_U_deformation(Tetrahedra& tet,CellCornerAttribute<vec3>& U) {
	
	// brush
	CellAttribute<int> dist(tet, std::numeric_limits<int>::max());
	for (auto seed : tet.iter_cells()) {
		if (dist[seed] < std::numeric_limits<int>::max()) continue;
		std::queue<int> file;
		dist[seed] = 0;
		file.push(seed);
		while (!file.empty()) {
			Volume::Cell cur(tet, file.front());
			file.pop();
			for (auto f : cur.iter_facets()) {
				auto opp = f.opposite();
				if (!opp.active()) continue;
				if (dist[opp.cell()] <= dist[cur] + 1) continue;
				GPTransitionFunction tf(f, U);
				for (auto corner : opp.cell().iter_corners()) U[corner] = tf.apply(U[corner]);
				dist[opp.cell()] = dist[cur] + 1;
				file.push(opp.cell());
			}
		}
	}

	// render
	Tetrahedra render_tet;
	render_tet.points.create_points(4 * tet.ncells());
	FOR(c, tet.ncorners()) render_tet.points[c] = U[c];
	render_tet.create_cells(tet.ncells());
	FOR(c, render_tet.ncorners()) render_tet.vert(c / 4, c % 4) = c;
	ToolBox(render_tet).drop_volume("deformed");
	//DropVolume(render_tet).apply("deformed");
}