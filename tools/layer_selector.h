#include <vector>
#include "../context.h"

struct LayerSelector {

	LayerSelector(Context &ctx) : ctx(ctx) {}

	void init();
	void hover();
	void select();
	void highlight_hovered_cells();
	void highlight_selected_cells();
	void clear();

	private:
	Context &ctx;
	std::vector<int> layers;

	int hovered_h = -1;

	std::vector<int> hovered_halfedges;
	// std::vector<int> last_hovered_halfedges;
	
	std::vector<int> selected_halfedges;
	// std::vector<int> last_selected_halfedges;

	int selected_layer = -1;

};