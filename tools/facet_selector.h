#include <vector>
#include "../context.h"

struct FacetSelector {

	FacetSelector(Context &ctx) : ctx(ctx) {}

	void paint(unsigned int hovered_facet);
	void clear();

	private:
	std::vector<int> last_hovered_facets;
	Context &ctx;
};