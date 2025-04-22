#include <vector>

struct SelectionCompressor {

	// Compress
	SelectionCompressor(std::vector<bool> &selection) : selection(selection) {
		for (int i = 0; i < selection.size(); ++i) {
			if (selection[i]) 
				compress_selection.push_back(i);
		}
	}

	// Decompress
	SelectionCompressor(std::vector<int> &compress_selection, int n) : compress_selection(compress_selection) {
		compress_selection.resize(n, false);
		for (auto x : compress_selection) {
			selection[x] = true;
		}
	}

	std::vector<int> get_compress_selection() {
		return compress_selection;
	}
	
	std::vector<bool> get_selection() {
		return selection;
	}

	private:
	std::vector<int> compress_selection;
	std::vector<bool> selection;

};