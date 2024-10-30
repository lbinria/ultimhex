#ifndef TOOLBOX_POINTSET_H__
#define TOOLBOX_POINTSET_H__


struct ToolBoxPointSet {
	ToolBoxPointSet(PointSet& pts);

	BBox3 bbox();

	BBox3 normalize(BBox3 bbox = BBox3());

	void normalize_for_stability();
	BBox3 set_min_to_zero();


	void merge_points(std::vector<int>& old2new, double epsilon = 0);

	void clear();
	void copy_from(PointSet& other);
	void backup(PointAttribute<vec3>& save);
	void restore(PointAttribute<vec3>& save);

	PointSet& pts;
};


template<>
struct ToolBox<PointSet> : public ToolBoxPointSet {
	ToolBox(PointSet& pts) :ToolBoxPointSet(pts) {  }
};
#endif