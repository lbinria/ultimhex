#include "mesh_metadata.h"

using json = nlohmann::json;

json MeshMetadata::MetadataAttribute::to_json() {

	json j;
	j["name"] = name;
	j["type"] = type;
	j["where"] = where;

	return j;
}

MeshMetadata::MetadataAttribute MeshMetadata::MetadataAttribute::from_json(json &j) {

	MetadataAttribute obj;
	obj.name = j["name"].get<std::string>();
	obj.type = j["type"].get<std::string>();
	obj.where = j["where"].get<GEO::MeshElementsFlags>();

	return obj;
}

json MeshMetadata::to_json() {

	std::vector<json> json_attributes(attributes.size());
	std::transform(attributes.begin(), attributes.end(), json_attributes.begin(), [](auto &x) {
		return x.to_json();
	});
	
	json j;
	j["filename"] = filename;
	j["tet_filename"] = tet_filename;
	j["cell_type"] = cell_type;
	j["attributes"] = json_attributes;
	
	return j;
}

MeshMetadata MeshMetadata::from_json(json &j) {

	MeshMetadata obj;
	obj.filename = j["filename"].get<std::string>();
	obj.tet_filename = j["tet_filename"].get<std::string>();
	obj.cell_type = j["cell_type"].get<GEO::MeshCellType>();

	std::vector<json> json_attributes = j["attributes"].get<std::vector<json>>();
	std::vector<MetadataAttribute> metadata_attributes(json_attributes.size());

	std::transform(json_attributes.begin(), json_attributes.end(), metadata_attributes.begin(), [](auto &sj) {
		return MeshMetadata::MetadataAttribute::from_json(sj);
	});

	obj.attributes = metadata_attributes;

	return obj;
}

bool MeshMetadata::has_attr(std::string attr_name) {
	for (auto attr : attributes) {
		if (attr.name == attr_name)
			return true;
	}

	return false;
}

std::optional<MeshMetadata::MetadataAttribute> MeshMetadata::get_attr(std::string attr_name) {
	for (auto attr : attributes) {
		if (attr.name == attr_name)
			return attr;
	}

	return std::nullopt;
}

void MeshMetadata::save() {
	// TODO check with filesystem::path if its valid path
	assert(filename != "");

	std::string j = to_json().dump();
	std::ofstream ofs(filename + ".json");
	ofs << j;
	ofs.close();
}