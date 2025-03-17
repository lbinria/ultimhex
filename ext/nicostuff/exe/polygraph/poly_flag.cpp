#include <ultimaille/all.h>
#include <framework/trace.h>

#include <framework/nico_framework.h>


using namespace UM::Linear;


// convention for flag direction -x,x,-y,y-z,z
std::array<vec3, 6> axis = { 
	vec3(-1,0,0), vec3(1,0,0), 
	vec3(0,-1,0), vec3(0,1,0), 
	vec3(0,0,-1), vec3(0,0,1) 
};


const int flag_step = 5;
const int flag_constraint_shift = 2;



void lockbend(Triangles& tri, FacetAttribute<int>& flag) {
	for (auto f : tri.iter_facets())
		flag[f] = flag_step * (flag[f] / flag_step) + flag_constraint_shift;
}
void unlockbend(Triangles& tri, FacetAttribute<int>& flag) {
	for (auto f : tri.iter_facets())
		flag[f] = flag_step * (flag[f] / flag_step) ;
}

void autofill(Triangles& tri, FacetAttribute<int>& flag) {
	auto flag_is_constained = [&](Surface::Facet f) {
		um_assert(flag[f] >= 0);
		return flag_constraint_shift  == flag[f] %flag_step;
	};

	auto flag_dim = [&](Surface::Facet f) {
		um_assert(flag[f] >= 0);
		int i = flag[f] / flag_step;
		return i / 2;
		};
	auto flag_sign = [&](Surface::Facet f) {
		um_assert(flag[f] >= 0);
		int i = flag[f] / flag_step;
		return (i % 2)  ? -1. : 1.;
	};
	FacetAttribute<mat3x3> B(tri, mat3x3::identity());
	mat3x3 I3 = mat3x3::identity();
	{
		LeastSquares ls(3 * tri.nfacets());
		bool have_constraint = false;
		for (auto f : tri.iter_facets()) {
			if (flag[f] < 0) continue;
			if (!flag_is_constained(f)) continue;
			if (Triangle3(f).unsigned_area() < 1e-10) continue;
			vec3 n = Triangle3(f).normal();

			vec3 rot = -cross(n, flag_sign(f)*I3[flag_dim(f)]);
			if (rot.norm2()>1e-15)
				rot = std::asin(std::min(1., rot.norm())) * rot.normalized();
			FOR(d, 3) ls.fix(3 * f + d, rot[d]);

		}
		for (auto h : tri.iter_halfedges()) {
			auto opp = h.opposite();
			if (!opp.active()) continue;
			FOR(d, 3) ls.add_to_energy(X(d + 3 * h.facet()) - X(d + 3 * opp.facet()));
		}
		ls.solve();
		for (auto f : tri.iter_facets()) {
			Quaternion q;
			FOR(d, 3) q.v[d] = ls.value(d + 3 * f);
			//plop(q.v);
			double alpha = q.v.norm();
			if (alpha < 1e-10) continue;// keep identity
			q.v = std::sin(alpha * .5) * q.v.normalized();
			q.w = std::cos(alpha * .5);
			B[f] = q.rotation_matrix().transpose();
		}
	}

	FacetAttribute<vec3> vd(tri);
	FOR(dim, 3) {
		for (auto f : tri.iter_facets()) vd[f] = B[f][dim];
		Drop(tri, vd).apply_line("vec");
	}

	for (auto f : tri.iter_facets()) {
		if (flag[f] >= 0 && flag_is_constained(f)) continue;
		vec3 n = Triangle3(f).normal();
		if ((n.norm2() - 1) > .1) n = vec3(1, 0, 0);
		int j = 0;
		FOR(i, 6) {
			double sign_i = (i % 2) ? (-1) : 1;
			int d_i = i / 2;
			double sign_j = (j % 2) ? (-1) : 1;
			int d_j = j / 2;
			if (sign_i * n * B[f][d_i] > sign_j * n * B[f][d_j])j = i;
		}
		flag[f] = flag_step * j;
	}
}




	void flag_diagnostic(Triangles & tri, FacetAttribute<int>&flag, std::string path) {
		auto flag_dim = [&](Surface::Facet f) {
		int i = flag[f] / flag_step;
		return i / 2;
		};
	auto flag_sign = [&](Surface::Facet f) {
		um_assert(flag[f] >= 0);
		int i = flag[f] / flag_step;
		return (i % 2) ? -1. : 1.;
	};

	for (auto f : tri.iter_facets()) um_assert(flag[f] >= 0);

	Trace::step("poisson like constrained opti");
	PointAttribute<vec3> pos(tri);
	FOR(dim, 3) {
		ConstrainedLeastSquares ls(tri.nverts());
		ls.add_to_constraints(X(0));
		for (auto f : tri.iter_facets()) 
			if (flag_dim(f) == dim) 
				for (auto h : f.iter_halfedges())
					ls.add_to_constraints(X(h.from()) - X(h.to()));	
		
		for (auto h : tri.iter_halfedges()) 
			ls.add_to_energy(X( h.from()) - X(h.to()) + (h.to().pos() - h.from().pos())[dim]);
		ls.solve();

		for (auto v : tri.iter_vertices()) pos[v][dim] = ls.value(v);
	}

	Trace::step("compute charts");
	FacetAttribute<int> chart(tri);
	{
		DisjointSet ds(tri.nfacets());
		for (auto h : tri.iter_halfedges()) {
			auto opp = h.opposite();
			um_assert(opp.active());
			if(flag[h.facet()] / flag_step == flag[opp.facet()] / flag_step )
				ds.merge(opp.facet(), h.facet());
		}
		ds.get_sets_id(chart.ptr->data);
	}
	Drop(tri, chart).apply("charts");

	Trace::step("find null frontiers");
	PolyLine pl;
	pl.points = tri.points;
	EdgeAttribute<int> fail_id(pl);
	int nfails = 0;
	auto segmentation_node = [&](Surface::Vertex v) {
		std::set<int> neig_charts;
		for (auto cir : v.iter_halfedges()) neig_charts.insert(chart[cir.facet()]);
		return neig_charts.size() > 2;
	};


	for (auto seed : tri.iter_halfedges()) {
		if (!segmentation_node(seed.from())) continue;
		if (chart[seed.facet()] == chart[seed.opposite().facet()]) continue;
		Surface::Halfedge cur = seed;
		std::vector<int> cur_ids;
		while (!segmentation_node(cur.to())) {
			cur_ids.push_back(cur);
			cur = cur.next(); 
			auto opp = cur.opposite();
			um_assert(opp.active());
			while (chart[cur.facet()] == chart[opp.facet()]) {
				cur = opp.next();
				opp = cur.opposite();
			}
		}
		cur_ids.push_back(cur);
		if ((pos[cur.to()] - pos[seed.from()]).norm2() < 1e-20) {
			for (int hid : cur_ids) {
				Surface::Halfedge h(tri, hid);
				int e = pl.create_edges(1);
				pl.vert(e, 0) = h.from();
				pl.vert(e, 1) = h.to();
				fail_id[e] = nfails;
			}
			nfails++;
		}
	}

	{// local issues
		// opposite flag touch each other
		for (auto h : tri.iter_halfedges()) {
			auto opp = h.opposite();
			um_assert(opp.active());
			if (flag_dim(h.facet()) != flag_dim(opp.facet())) continue;
			if (flag_sign(h.facet()) != flag_sign(opp.facet())) {
				int e = pl.create_edges(1);
				pl.vert(e, 0) = h.from();
				pl.vert(e, 1) = h.to();
				fail_id[e] = -2;
			}
		}
	}

	Drop(pl, fail_id)._skip_value(-10).apply("fail");
	DropPolyLine(pl).add(fail_id,"fail_id")._just_save_filename(path + "/feedback.geogram").apply();
}

void framework_parameters(NicoFramework& fw) {
	fw.add("enum", "ALGO", "lockbend").possible_values("create,fill,diagnostic,lockbend,unlockbend").description("algorithm to run");
	//fw.add("string", "projectpath", "C:/NICO/data/polycubeprojects/pain");
	fw.add("string", "projectpath", "C:/NICO/data/polycubeprojects/CAD/B21");
	//fw.add("string", "projectpath", "C:/NICO/data/polycubeprojects/easy_pont");
	
}
void framework_main(NicoFramework& fw) {
	plop(std::string(fw["projectpath"]));
	
	std::string path = fw["projectpath"];
	Triangles tri;
	auto attr = read_by_extension(path + "/flag.geogram", tri);

	if (tri.nfacets() == 0) {
		Trace::step("Flag field does not exist, create a new one");
		Triangles m;
		read_by_extension(path + "/tri.geogram", m);
		{
			// try to fix common incompatibilities
			ToolBox(m).make_orientable_manifold();
			double vol = 0;
			//m.connect();
			for (auto f : m.iter_facets()) vol += Tetrahedron(vec3(0, 0, 0), f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos()).volume();
			plop(vol);
			if (vol < 0) {
				m.disconnect();
				for (auto f : m.iter_facets())
					std::swap(m.vert(f, 1), m.vert(f, 2));
			}
		}
		FacetAttribute<int> flag(m, -1);

		//DropSurface(m).add(flag, "flag").apply();return;
		
		DropSurface(m).add(flag, "flag")._just_save_filename(path + "/flag.geogram").apply();
		attr = read_by_extension(path + "/flag.geogram", tri);
	}
	tri.connect();
	FacetAttribute<int> flag("flag", attr, tri);
	
	if (fw["ALGO"].is("diagnostic")) {
		flag_diagnostic(tri, flag, path);
		return;
	}
	
	


	Drop(tri, flag).apply("inputflag");
	if (fw["ALGO"].is("unlockbend"))  unlockbend(tri, flag);
	if (fw["ALGO"].is("lockbend"))  lockbend(tri, flag);
	if (fw["ALGO"].is("fill"))      autofill(tri, flag);
	Drop(tri, flag).apply("outputflag");

	if (fw["run_from"].is("graphite")) 
	DropSurface(tri).add(flag, "flag")._just_save_filename(path + "/flag.geogram").apply();

}
