#include <toolbox.h>

	ToolBoxPointSet::ToolBoxPointSet(PointSet& pts) :pts(pts) {  }

	BBox3 ToolBoxPointSet::bbox() {
		BBox3 res;
		FOR(i, pts.size()) res.add(pts[i]);
		return res;
	}

	
	
	void ToolBoxPointSet::normalize_largest_dimension() {
		auto b = bbox();
		double inv_factor = 1e20;
		FOR(d, 3) inv_factor = std::min(inv_factor, b.size()[d] );
		FOR(d, 3) FOR(i, pts.size()) pts[i][d] = inv_factor *(pts[i][d]-b.min[d]);
		
	}

	BBox3 ToolBoxPointSet::normalize(BBox3 bbox ) {
		if (bbox.empty()) FOR(i, 2) bbox.add(vec3(i, i, i));
		BBox3 prev_bbox;
		FOR(i, pts.size()) prev_bbox.add(pts[i]);
		double factor = 1e20;
		FOR(d, 3) {
			factor = (bbox.max[d] - bbox.min[d]) / (prev_bbox.max[d] - prev_bbox.min[d]);
			FOR(i, pts.size()) pts[i][d] = bbox.min[d] + factor * (pts[i][d] - prev_bbox.min[d]);
		}
		return prev_bbox;
	}

	void ToolBoxPointSet::normalize_for_stability() {
		BBox3 bbox;
		FOR(i, 2) bbox.add(10 * vec3(-1, -1, -1) + 20 * vec3(i, i, i));

		BBox3 prev_bbox;
		FOR(i, pts.size()) prev_bbox.add(pts[i]);
		double factor = 1e20;
		FOR(d, 3) factor = std::min(factor, (bbox.max[d] - bbox.min[d]) / (prev_bbox.max[d] - prev_bbox.min[d]));

		FOR(d, 3) FOR(i, pts.size()) pts[i][d] = bbox.min[d] + factor * (pts[i][d] - prev_bbox.min[d]);
	}
	BBox3 ToolBoxPointSet::set_min_to_zero() {
		BBox3 prev_bbox;
		FOR(i, pts.size()) prev_bbox.add(pts[i]);
		FOR(i, pts.size()) pts[i] = (pts[i] - prev_bbox.min);
		return prev_bbox;
	}


	void ToolBoxPointSet::merge_points(std::vector<int>& old2new, double epsilon ) {

		UM::colocate(*(pts.data), old2new, epsilon);
		std::vector<bool> tokill(pts.size());
		FOR(v, pts.size()) tokill[v] = (old2new[v] != v);
		std::vector<int> old2new2;
		pts.delete_points(tokill, old2new2);
		for (int& it : old2new) it = old2new2[it];
	}

	void ToolBoxPointSet::clear() { std::vector<bool> to_kill(pts.size(), true); std::vector<int> trashmap; pts.delete_points(to_kill, trashmap); }
	void ToolBoxPointSet::copy_from(PointSet& other) { if (pts.size() != 0) clear();  pts.create_points(other.size()); FOR(v, pts.size()) pts[v] = other[v]; }
	void ToolBoxPointSet::backup(PointAttribute<vec3>& save) { FOR(v, pts.size()) save[v] = pts[v]; }
	void ToolBoxPointSet::restore(PointAttribute<vec3>& save) { FOR(v, pts.size()) pts[v] = save[v]; }



