#include <ultimaille/all.h>
#include <create_shape.h>
#include <basic.h>
#include <drop_attribute.h>
#include <toolbox_hexahedra.h>
using namespace UM;


void CreatePointSet::from_vector(PointSet& p, const std::vector<vec3>& data) {
	p.create_points(data.size());
	FOR(i, data.size()) p[i] = data[i];
}

void CreatePointSet::grid1d(PointSet& p, int n) {
	p.create_points(n);
	FOR(i, n) p[i] = (1. / double(n)) * vec3(i, 0, 0);
}

void CreatePointSet::grid2d(PointSet& p,int n) {
	p.create_points(n * n);
	FOR(i, n * n) p[i] = (1. / double(n)) * vec3(i / n, i % n, 0);
}
void CreatePointSet::grid3d(PointSet& p,int n) {
	p.create_points(n * n * n);
	int n2 = n * n;
	FOR(i, n * n * n) p[i] = (1. / double(n)) * vec3(i / n2, (i - (i / n2) * n2) / n, i % n);
}



void CreatePolyline::grid1d(PolyLine& p,int n) {
	CreatePointSet::grid1d(p.points ,n);
	p.create_edges(n - 1);
	FOR(i, n - 1) {
		p.vert(i, 0) = i;
		p.vert(i, 1) = i + 1;
	}
}
void CreatePolyline::grid2d(PolyLine& p, int n) {
	CreatePointSet::grid2d(p.points,n);
	p.create_edges(2 * n * (n - 1));
	int seg = 0;
	FOR(i, n - 1) FOR(j, n) {
		p.vert(seg, 0) = i + j * n;
		p.vert(seg, 1) = i + 1 + j * n;
		seg++;
	}
	FOR(i, n) FOR(j, n - 1) {
		p.vert(seg, 0) = j * n + i;
		p.vert(seg, 1) = (j + 1) * n + i;
		seg++;
	}
}
void CreatePolyline::circle(PolyLine& p, int n) {
	p.points.create_points(n);
	FOR(v, n) p.points[v] = vec3(cos(2. * M_PI * double(v) / double(n)), sin(2. * M_PI * double(v) / double(n)), 0);
	p.create_edges(n);
	FOR(i, n) {
		p.vert(i, 0) = i;
		p.vert(i, 1) = (i + 1) % n;
	}
}


template<class MESH>
void set_vertices(MESH& m,int element, const std::vector<int>& data) {
	FOR(i, data.size()) m.vert(element, i) = data[i];
}

void CreateSurface::quad_grid(Quads& m,int n) {
	CreatePointSet::grid2d(m.points,n + 1);
	m.create_facets(n * n);
	int f = 0;
	FOR(i, n) FOR(j, n) {
		set_vertices(m, f, {
			(i + 0) + (j + 0) * (n + 1),
			(i + 1) + (j + 0) * (n + 1),
			(i + 1) + (j + 1) * (n + 1),
			(i + 0) + (j + 1) * (n + 1)
		});
		f++;
	}
}
void CreateSurface::triangles_grid(Triangles& m, int n) {
	CreatePointSet::grid2d(m.points,n + 1);
	m.create_facets(2 * n * n);
	int f = 0;
	FOR(i, n) FOR(j, n) {
		m.vert(f, 0) = (i + 0) + (j + 0) * (n + 1);
		m.vert(f, 1) = (i + 1) + (j + 0) * (n + 1);
		m.vert(f, 2) = (i + 1) + (j + 1) * (n + 1);
		f++;
		m.vert(f, 0) = (i + 0) + (j + 0) * (n + 1);
		m.vert(f, 1) = (i + 1) + (j + 1) * (n + 1);
		m.vert(f, 2) = (i + 0) + (j + 1) * (n + 1);
		f++;
	}
}

void CreateSurface::tet(Triangles& m) {
	CreatePointSet::from_vector(m.points, { vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1)});
	m.create_facets(4);
	set_vertices(m, 0, { 0,1,2 });
	set_vertices(m, 1, { 1,0,3 });
	set_vertices(m, 2, { 2,1,3 });
	set_vertices(m, 3, { 3,0,2 });
}


void CreateSurface::non_manifold_vertex1(Triangles& m) {
	m.points.create_points(13);
	FOR(v, 6) m.points[v] = vec3(cos(2. * M_PI * double(v) / 6.), sin(2. * M_PI * double(v) / 6.), -.5);
	FOR(v, 6) m.points[v + 6] = m.points[v] + vec3(0, 0, 1);
	m.create_facets(12);
	FOR(f, 6) {
		set_vertices(m, f, { 12,f, (f + 1) % 6 });
		set_vertices(m, f+6, { 12,f+6, ((f + 1) % 6) + 6 });
	}
}
void CreateSurface::non_manifold_vertex2(Triangles& m) {
	CreatePointSet::from_vector(m.points, { vec3(0,0,0), vec3(1,0,0), vec3(1,1,0), vec3(-1,0,0), vec3(-1,1,0) });
	m.create_facets(2);
	set_vertices(m, 0, { 0,1,2 });
	set_vertices(m, 1, { 0,3,4 });
}
void CreateSurface::non_manifold_edge(Triangles& m) {
	CreatePointSet::from_vector(m.points, { vec3(0,0,0), vec3(1,0,0), vec3(1,1,0), vec3(1,1,1)});
	m.create_facets(2);
	set_vertices(m, 0, { 0,1,2 });
	set_vertices(m, 1, { 0,1,3 });
}


void CreateSurface::triangles_fan_negative_curvature(Triangles& m) {
	m.points.create_points(13);
	FOR(v, 12) m.points[v] = vec3(cos(2. * M_PI * double(v) / 12.), sin(2. * M_PI * double(v) / 12.), .5*double(-1+2*(v%2)));
	m.create_facets(12);
	FOR(f, 12) set_vertices(m, f, { 12,f,(f + 1) % 12 });
}


void CreateSurface::polygonal_hut(Polygons& m, double radius, double height, int n) {
	m.points.create_points(2 * n + 1);
	FOR(v, n) m.points[v] = radius * vec3(cos(2. * M_PI * double(v) / double(n)), sin(2. * M_PI * double(v) / double(n)), 0);
	FOR(v, n) m.points[v + n] = m.points[v] + vec3(0, 0, height);
	m.points[2 * n] = vec3(0, 0, height + radius);
	int f = m.create_facets(n, 4);
	FOR(i, n) {
		m.vert(f, 0) = (i + 0) % n + 0;
		m.vert(f, 1) = (i + 1) % n + 0;
		m.vert(f, 2) = (i + 1) % n + n;
		m.vert(f, 3) = (i + 0) % n + n;
		f++;
	}
	f = m.create_facets(n, 3);
	FOR(i, n) {
		m.vert(f, 0) = (i + 0) % n + n;
		m.vert(f, 1) = (i + 1) % n + n;
		m.vert(f, 2) = 2 * n;
		f++;
	}
}






void CreateVolume::edge_one_ring(Tetrahedra& m,int n ) {
	m.points.create_points(n + 2);
	FOR(v, n) m.points[v] = vec3(cos(2. * M_PI * double(v) / double(n)), sin(2. * M_PI * double(v) / double(n)), 0);
	m.points[n] = vec3(0, 0, 1);
	m.points[n + 1] = vec3(0, 0, -1);
	m.create_cells(n);
	FOR(i, n) {
		m.vert(i, 0) = i;
		m.vert(i, 1) = (i + 1) % n;
		m.vert(i, 2) = n;
		m.vert(i, 3) = n + 1;
	}
}

void CreateVolume::hex_grid(Hexahedra&m, int n) {
	int n1 = n + 1;
	CreatePointSet::grid3d(m.points,n1);
	m.create_cells(n * n * n);
	FOR(i, n)FOR(j, n)FOR(k, n) {
		int c = i + n * j + n * n * k;
		FOR(di, 2)FOR(dj, 2)FOR(dk, 2)
			m.vert(c, di + 2 * dj + 4 * dk) = i + di + n1 * (j + dj) + n1 * n1 * (k + dk);
	}
}

void CreateVolume::tri_grid(Hexahedra& grid, int nb_subdivision) {
	grid.points.create_points(14);
	grid.points[0] = vec3(0, 0, 0);
	grid.points[1] = vec3(2, 0, 0);
	grid.points[2] = vec3(3, 2, 0);
	grid.points[3] = vec3(2, 4, 0);
	grid.points[4] = vec3(0, 4, 0);
	grid.points[5] = vec3(-2, 2, 0);
	grid.points[6] = vec3(1, 2, 0);
	FOR(i, 7) grid.points[7 + i] = grid.points[i] + vec3(0, 0, 2);
	grid.create_cells(3);
	int cell_bottom_vertices[3][4] = { {0,1,6,2},{2,3,6,4},{4,5,6,0} };
	FOR(c, 3) FOR(lv, 4) grid.vert(c, lv) = cell_bottom_vertices[c][lv];
	FOR(c, 3) FOR(lv, 4) grid.vert(c, lv + 4) = cell_bottom_vertices[c][lv] + 7;
	FOR(i, nb_subdivision) ToolBox(grid).split();
}


void CreateVolume::two_wedges(Wedges& m) {
	CreatePointSet::grid3d(m.points,2);
	m.create_cells(2);
	int data[2][6] = { {0,1,2,4,5,6},{3,2,1,7,6,5} };
	FOR(c, 2)FOR(lv, 6) m.vert(c, lv)=data[c][lv];
}

void CreateVolume::pyramids_cube(Pyramids& m) {
	m.points.create_points(9) ;
	FOR(i, 2)FOR(j, 2)FOR(k, 2) m.points[i+2*j+4*k] = vec3(i ,j,k);
	m.points[8] = vec3(.5, .5, .5);
	m.create_cells(6);
	int data[6][5] = { 
			{0,1,3,2,8},
			{1,5,7,3,8},
			{0,4,5,1,8},
			{4,6,7,5,8},
			{2,3,7,6,8},
			{0,2,6,4,8}
	};
	FOR(c, 6)FOR(lv, 5) m.vert(c, lv) = data[c][lv];
}



void test_create_shapes() {
	PointSet p0;
	CreatePointSet::grid1d(p0, 5);
	DropPointSet(p0).apply("point_grid_1d");

	PointSet p1;
	CreatePointSet::grid2d(p1, 5);
	DropPointSet(p1).apply("point_grid_2d");

	PointSet p2;
	CreatePointSet::grid3d(p2, 5);
	DropPointSet(p2).apply("point_grid_3d");


	PolyLine s0;
	CreatePolyline::grid1d(s0, 5);
	DropPolyLine(s0).apply("polyline_grid_1d");

	PolyLine s1;
	CreatePolyline::grid2d(s1, 5);
	DropPolyLine(s1).apply("polyline_grid_1d");

	PolyLine s2;
	CreatePolyline::circle(s2, 20);
	DropPolyLine(s2).apply("polyline_circle");


	Triangles t0;
	CreateSurface::triangles_grid(t0, 10);
	DropSurface(t0).apply("triangles_grid_2d");

	Triangles t1;
	CreateSurface::triangles_fan_negative_curvature(t1);
	DropSurface(t1).apply("triangles_fan");

	Triangles t2;
	CreateSurface::non_manifold_vertex1(t2);
	DropSurface(t2).apply("non_manifold_vertex1");

	Triangles t3;
	CreateSurface::non_manifold_vertex2(t3);
	DropSurface(t3).apply("non_manifold_vertex2");

	Triangles t4;
	CreateSurface::non_manifold_edge(t4);
	DropSurface(t4).apply("non_manifold_edge");

	Triangles t5;
	CreateSurface::tet(t5);
	DropSurface(t5).apply( "tet");





	Quads q;
	CreateSurface::quad_grid(q, 10);
	DropSurface(q).apply("quads_grid_2d");
	
	Polygons po;
	CreateSurface::polygonal_hut(po);
	DropSurface(po).apply("hut_2d");

	Tetrahedra te;
	CreateVolume::edge_one_ring(te);
	DropVolume(te).apply("tet_one_ring");
	
	Hexahedra h;
	CreateVolume::hex_grid(h, 3);
	DropVolume(h).apply("hex_grid");

	Wedges w;
	CreateVolume::two_wedges(w);
	DropVolume(w).apply("wedge_grid");

	Pyramids py;
	CreateVolume::pyramids_cube(py);
	DropVolume(py).apply("pyra_grid");


}