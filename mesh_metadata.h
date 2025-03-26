#pragma once

#include <optional>

#include <geogram/mesh/mesh.h>
#include <json.hpp>

using json = nlohmann::json;

struct MeshMetadata {

	struct MetadataAttribute {
		std::string name;
		std::string type;
		GEO::MeshElementsFlags where;
		// std::optional<std::pair<double, double>> range;
		json to_json();
		static MetadataAttribute from_json(json &json);

	};

	std::string filename;
	std::string tet_filename;
	GEO::MeshCellType cell_type;
	std::vector<MetadataAttribute> attributes;

	json to_json();
	static MeshMetadata from_json(json &json);
	bool has_attr(std::string attr_name);
	std::optional<MetadataAttribute> get_attr(std::string attr_name);
	void save();

};