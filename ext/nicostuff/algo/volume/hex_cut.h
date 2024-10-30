#ifndef HEX_CUT__H__
#define	HEX_CUT__H__

#include <ultimaille/all.h>
#include <iostream>
#include <algo/volume/SDF.h>

struct HexCutter{
	HexCutter(Hexahedra& hex);
	~HexCutter();

	void init_from_charts();

	void add_sdf(SDF* s);

	void show();

	void show_feature_neigs(std::string name = "featureneig", bool include_boundary = false, bool only_last_SDF = false);


	Hexahedra& hex;
	//EC3d conn;
	CellFacetAttribute<int> chart;
	std::vector<SDF*> sdf;
	bool drop_all_steps;
private:
	void optimize_SDF(int sdf_id);
};
#endif
