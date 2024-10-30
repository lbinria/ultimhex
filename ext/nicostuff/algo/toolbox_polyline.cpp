#include<algo/toolbox.h>



ToolBoxPolyLine::ToolBoxPolyLine(PolyLine& pl) :pl(pl) {  }

	void ToolBoxPolyLine::triangulate_with_shewchuck(Triangles& m, int verbose , int may_add_vertices ) {

		struct triangulateio in, out;

		in.numberofpoints = pl.nverts();
		in.pointlist = (REAL*)malloc(in.numberofpoints * 2 * sizeof(REAL));
		FOR(v, in.numberofpoints) FOR(d, 2) in.pointlist[2 * v + d] = pl.points[v][d];

		in.numberofsegments = pl.nedges();
		in.segmentlist = (int*)malloc(2 * in.numberofsegments * sizeof(int));
		FOR(seg, in.numberofsegments) FOR(e, 2) in.segmentlist[2 * seg + e] = pl.vert(seg, e);


		in.numberofpointattributes = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.numberoftriangles = 0;


		in.pointmarkerlist = (int*)NULL; /* Not needed if -N or -B switch used. */
		in.segmentmarkerlist = (int*)NULL;

		out.pointlist = (REAL*)NULL;
		out.trianglelist = (int*)NULL;
		out.pointattributelist = (REAL*)NULL;
		out.pointmarkerlist = (int*)NULL; /* Not needed if -N or -B switch used. */
		/* Not needed if -E switch used or number of triangle attributes is zero: */
		out.triangleattributelist = (REAL*)NULL;
		out.neighborlist = (int*)NULL;         /* Needed only if -n switch used. */
		/* Needed only if segments are output (-p or -c) and -P not used: */
		out.segmentlist = (int*)NULL;
		/* Needed only if segments are output (-p or -c) and -P and -B not used: */
		out.segmentmarkerlist = (int*)NULL;
		out.edgelist = (int*)NULL;             /* Needed only if -e switch used. */
		out.edgemarkerlist = (int*)NULL;   /* Needed if -e used and -B not used. */

		plop(verbose);
		plop(may_add_vertices);
		if (may_add_vertices) {

			double a = ave_edge_size();
			a *= a;
			char data[1024];
			sprintf(data, "pa%fCVz", a);
			if (verbose)	triangulate(data, &in, &out, NULL);
			else			triangulate((char*)("paQz"), &in, &out, NULL);
		}
		else {
			if (verbose)	triangulate((char*)("pCVz"), &in, &out, NULL);
			else			triangulate((char*)("pQz"), &in, &out, NULL);
		}//if (out.numberoftriangles > 100000) return;

		m.points.create_points(out.numberofpoints);
		FOR(v, out.numberofpoints) FOR(d, 2) m.points[v][d] = out.pointlist[v * 2 + d];

		m.create_facets(out.numberoftriangles);
		//plop(out.numberoftriangles);
		FOR(f, m.nfacets()) FOR(lv, 3) m.vert(f, lv) = out.trianglelist[3 * f + lv];

		free(in.pointlist);
		free(in.segmentlist);
		free(out.pointlist);
		free(out.trianglelist);

	}
	double ToolBoxPolyLine::ave_edge_size() {
		double sum = 0;
		FOR(s, pl.nedges()) sum += (pl.points[pl.vert(s, 0)] - pl.points[pl.vert(s, 1)]).norm();
		return sum / double(pl.nedges());
	}
	void ToolBoxPolyLine::kill_isolated_vertices() {
		std::vector<bool> to_kill(pl.nverts(), true);
		FOR(s, pl.nedges()) FOR(lv, 2) to_kill[pl.vert(s, lv)] = false;
		pl.delete_vertices(to_kill);
	}
	int ToolBoxPolyLine::add_segment(vec3 A, vec3 B) {
		int s = pl.create_edges(1);
		int v = pl.points.create_points(2);
		FOR(lv, 2) pl.vert(s, lv) = v + lv;
		pl.points[v] = A;
		pl.points[v + 1] = B;
		return s;
	}
