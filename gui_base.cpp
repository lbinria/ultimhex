#include <geogram/basic/command_line.h> // for CmdLine::get_arg() and CmdLine::set_arg()

#include "gui_base.h"

SimpleMeshApplicationExt::ColorArray::ColorArray(std::initializer_list<std::array<float,4>> init_values) {
    chars_.resize(init_values.size()*4);
    floats_.resize(init_values.size()*4);
    std::size_t color_index = 0;
    for(auto &color : init_values) {
        floats_[4*color_index+0] = color[0]; // R
        floats_[4*color_index+1] = color[1]; // G
        floats_[4*color_index+2] = color[2]; // B
        floats_[4*color_index+3] = color[3]; // A
        update_chars_of_color(color_index);
        color_index++;
    }
}

float* SimpleMeshApplicationExt::ColorArray::as_floats() {
    return floats_.data();
}

float* SimpleMeshApplicationExt::ColorArray::color_as_floats(std::size_t color_index) {
    return (floats_.data()+(4*color_index));
}

ImVec4 SimpleMeshApplicationExt::ColorArray::color_as_ImVec4(std::size_t color_index) {
    return ImVec4(
        floats_.data()[4*color_index+0],
        floats_.data()[4*color_index+1],
        floats_.data()[4*color_index+2],
        floats_.data()[4*color_index+3]
    );
}

unsigned char* SimpleMeshApplicationExt::ColorArray::as_chars() {
    return chars_.data();
}

unsigned char* SimpleMeshApplicationExt::ColorArray::color_as_chars(std::size_t color_index) {
    return (chars_.data()+(4*color_index));
}

void SimpleMeshApplicationExt::ColorArray::update_chars_of_color(std::size_t color_index) {
    chars_[4*color_index+0] = static_cast<unsigned char>(floats_[4*color_index+0]*255.0f); // R
    chars_[4*color_index+1] = static_cast<unsigned char>(floats_[4*color_index+1]*255.0f); // G
    chars_[4*color_index+2] = static_cast<unsigned char>(floats_[4*color_index+2]*255.0f); // B
    chars_[4*color_index+3] = static_cast<unsigned char>(floats_[4*color_index+3]*255.0f); // A
}

SimpleMeshApplicationExt::SimpleMeshApplicationExt(const std::string &name) : SimpleMeshApplication(name) {}

SimpleMeshApplicationExt::~SimpleMeshApplicationExt() {
    clear_scene_overlay();
}


void SimpleMeshApplicationExt::geogram_initialize(int argc, char** argv) {
    SimpleMeshApplication::geogram_initialize(argc,argv);
#ifdef WINDOW_SIZE
    // Change window size
    // By order of priority:
    // - command line argument "gfx:geometry"
    // - preprocessor macro WINDOW_SIZE defined in CMakeLists.txt
    // - default value "1024x1024", hard-coded in ext/geogram/src/lib/geogram/basic/command_line_args.cpp import_arg_group_gfx()
    if(!phone_screen_ && CmdLine::get_arg("gfx:geometry") == "1024x1024") { // if not a phone screen and default value for gfx:geometry 

        CmdLine::set_arg("gfx:geometry", WINDOW_SIZE); // bigger window

    }
#endif
}

void SimpleMeshApplicationExt::clear_scene_overlay() {
    points_groups_.clear();
    edges_groups_.clear();
}

std::size_t SimpleMeshApplicationExt::new_points_group(const float* rgba, const float* size, bool* show) {
    points_groups_.push_back(PointsGroup(rgba,size,show));
    return index_of_last(points_groups_);
}

std::size_t SimpleMeshApplicationExt::new_edges_group(unsigned int colormap_index, double texture_coordinate, int* width, bool* show) {
    edges_groups_.push_back(EdgesGroup(colormap_index,texture_coordinate,width,show));
    return index_of_last(edges_groups_);
}

void SimpleMeshApplicationExt::add_point_to_group(std::size_t group, double x, double y, double z) {
    points_groups_.at(group).points.push_back({x,y,z});
}

void SimpleMeshApplicationExt::add_edge_to_group(std::size_t group, double x1, double y1, double z1, double x2, double y2, double z2) {
    edges_groups_.at(group).edges.push_back(
        std::make_pair<vec3,vec3>({x1,y1,z1},{x2,y2,z2})
    );
}

void SimpleMeshApplicationExt::set_points_group_color(std::size_t index, const float* new_color) {
    points_groups_.at(index).color = new_color;
}

void SimpleMeshApplicationExt::set_edges_group_color(std::size_t index, unsigned int new_colormap_index, double new_texture_coordinate) {
    edges_groups_.at(index).colormap_index = new_colormap_index;
    edges_groups_.at(index).texture_coordinate = new_texture_coordinate;
}

void SimpleMeshApplicationExt::draw_scene() {
    SimpleMeshApplication::draw_scene();

    // draw points groups
    for(const PointsGroup& group : points_groups_) {
        if(*group.show) {
            geo_assert(group.color != nullptr);
            glupSetColor4fv(GLUP_FRONT_COLOR, group.color); // use the color of this group of points
            glupSetPointSize(*group.size);
            for(std::size_t i = 0; i<group.points.size(); ++i) { // for each point in this group
                glupBegin(GLUP_POINTS);
                glupPrivateVertex3dv(group.points[i].data()); // give a pointer to the coordinates to GLUP
                glupEnd();
            }
        }
    }

    // draw edges groups
    for(const auto& group : edges_groups_) {
        if(*group.show) {
            glupEnable(GLUP_TEXTURING);
            glActiveTexture(GL_TEXTURE0 + GLUP_TEXTURE_2D_UNIT);
            glBindTexture(GL_TEXTURE_2D, (GLuint) colormaps_[group.colormap_index].texture);
            glupTextureType(GLUP_TEXTURE_2D);
            glupTextureMode(GLUP_TEXTURE_REPLACE);
            glupPrivateTexCoord1d(group.texture_coordinate);
            glupSetMeshWidth(*group.width);
            glupBegin(GLUP_LINES);
            for(const auto& edge : group.edges) { // for each edge in this group
                glupPrivateVertex3dv(edge.first.data()); // give a pointer to the coordinates to GLUP
                glupPrivateVertex3dv(edge.second.data()); // give a pointer to the coordinates to GLUP
            }
            glupDisable(GLUP_TEXTURING);
            glupEnd();
        }
    }
}

void SimpleMeshApplicationExt::cursor_pos_callback( double x, double y, int source ) {
    cursor_pos_ = GEO::vec2(x,y); // store locally the current cursor position
    SimpleMeshApplication::cursor_pos_callback(x,y,source);
}

index_t SimpleMeshApplicationExt::pick(MeshElementsFlags what) {

    // see https://github.com/BrunoLevy/geogram/discussions/88
    // see https://github.com/BrunoLevy/geogram/pull/102
    // based on https://github.com/BrunoLevy/GraphiteThree/blob/main/src/lib/OGF/renderer/context/rendering_context.cpp get_picked_point()

    index_t x = index_t(cursor_pos_.x), y = index_t(cursor_pos_.y); // double to integer conversion of current cursor position
    if(x >= get_width() || y >= get_height()) { // if cursor out of the window
        return index_t(-1);
    }
    y = get_height()-1-y; // change Y axis orientation. glReadPixels() wants pixel coordinates from bottom-left corner

    mesh_gfx()->set_picking_mode(what); // instead of rendering colors, mesh_gfx will render indices
    draw_scene(); // rendering
    // read the index of the picked element using glReadPixels()
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ROW_LENGTH, 1);
    glReadPixels(
        GLint(x),GLint(y),1,1,GL_RGBA,GL_UNSIGNED_BYTE,buffer
    );

	//  Read the z coordinate of the picked point
	// in the depth buffer.
	float z;
	glReadPixels(
		GLint(x), GLint(y),1,1, GL_DEPTH_COMPONENT, GL_FLOAT, &z
	);
	picked_depth_ = double(z);

	// Restore default value
	glPixelStorei(GL_PACK_ROW_LENGTH, 0);

    // Transform the picked point into world space.
    picked_point_ = unproject(vec3(x, y, picked_depth_));


    mesh_gfx()->set_picking_mode(MESH_NONE); // go back to color rendering mode
    // decode index from pixel color
    return index_t(buffer[0])        |
            (index_t(buffer[1]) << 8)  |
            (index_t(buffer[2]) << 16) |
            (index_t(buffer[3]) << 24);
}


std::vector<index_t> SimpleMeshApplicationExt::pick_size(MeshElementsFlags what, int size) {
    // Memory::byte buffer[size][size][4]; // a 4-bytes buffer to read pixels
    unsigned char buffer[size][size][4]; // a 4-bytes buffer to read pixels

    // see https://github.com/BrunoLevy/geogram/discussions/88
    // see https://github.com/BrunoLevy/geogram/pull/102
    // based on https://github.com/BrunoLevy/GraphiteThree/blob/main/src/lib/OGF/renderer/context/rendering_context.cpp get_picked_point()
    std::vector<index_t> picked_elements(size * size, -1);

    mesh_gfx()->set_picking_mode(what); // instead of rendering colors, mesh_gfx will render indices
    draw_scene(); // rendering
    // read the index of the picked element using glReadPixels()
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ROW_LENGTH, 1);

    index_t x = index_t(cursor_pos_.x - size / 2), y = index_t(cursor_pos_.y); // double to integer conversion of current cursor position
    if(x >= get_width() || y >= get_height()) { // if cursor out of the window
        return std::vector<index_t>();
    }
    y = get_height()-1-y; // change Y axis orientation. glReadPixels() wants pixel coordinates from bottom-left corner
    y -= size / 2;

    glReadPixels(
        GLint(x),GLint(y),size,size,GL_RGBA,GL_UNSIGNED_BYTE,buffer
    );
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            // decode index from pixel color
            index_t picked_element = index_t(buffer[j][i][0])        |
                    (index_t(buffer[j][i][1]) << 8)  |
                    (index_t(buffer[j][i][2]) << 16) |
                    (index_t(buffer[j][i][3]) << 24);

            picked_elements[i + j * size] = picked_element;
        }
    }


	// Restore default value
	glPixelStorei(GL_PACK_ROW_LENGTH, 0);

    mesh_gfx()->set_picking_mode(MESH_NONE); // go back to color rendering mode
    
    return picked_elements;
}

index_t SimpleMeshApplicationExt::pickup_cell_edge(GEO::vec3 p0, index_t c_idx) {
	// Search nearest edge
	double min_dist = std::numeric_limits<double>().max();
	index_t e_idx = -1;

	// Search nearest edge
	for (int i = 0; i < mesh_.cells.nb_edges(c_idx); i++) {
		// Get points from current edge
		index_t v0_idx = mesh_.cells.edge_vertex(c_idx, i, 0);
		index_t v1_idx = mesh_.cells.edge_vertex(c_idx, i, 1);
		GEO::vec3 &p1 = mesh_.vertices.point(v0_idx);
		GEO::vec3 &p2 = mesh_.vertices.point(v1_idx);

		// Compute dist from mouse point to edge points
		double dist = distance(p0, (p1 + p2) / 2.);

		// Keep min dist
		if (dist < min_dist) {
			min_dist = dist;
			e_idx = i;
			// posA = p1;
			// posB = p2;
		}
	}

	return e_idx;
}

std::tuple<index_t, index_t> SimpleMeshApplicationExt::pickup_cell_facet2(GEO::vec3 p0, index_t c_idx) {
	// Search if point is on facet
	double min_dist = std::numeric_limits<double>().max();
	index_t f_idx = NO_FACET;
	index_t lf_idx = NO_FACET;

	for (index_t lf = 0; lf < mesh_.cells.nb_facets(c_idx); lf++) {

		auto a = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 0));
		auto b = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 1));
		auto c = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 2));

		auto bary = (a + b + c) / 3.;

		auto n = normalize(cross(b - a, c - b));
		double dist = dot(p0 - a, n);
		
		if (std::abs(dot(normalize(p0 - a), normalize(bary - a))) < 1e-4)
			continue;

		if (dist < min_dist) {
			min_dist = dist;
			f_idx = mesh_.cells.facet(c_idx, lf);
			lf_idx = lf;
		}
	}

	return std::tuple<index_t, index_t>(f_idx, lf_idx);
}

// TODO improve by breaking loop when dist is less than a given threshold
std::tuple<index_t, index_t> SimpleMeshApplicationExt::pickup_cell_facet(GEO::vec3 p0, index_t c_idx) {
	// Search if point is on facet
	double min_dist = std::numeric_limits<double>().max();
	index_t f_idx = NO_FACET;
	index_t lf_idx = NO_FACET;

    
	for (index_t lf = 0; lf < mesh_.cells.nb_facets(c_idx); lf++) {

        int nb_verts = mesh_.cells.facet_nb_vertices(c_idx, lf);

		auto a = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 0));
		auto b = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 1));
		auto c = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 2));
        
        auto bary = a + b + c;

        for (int i = 3; i < nb_verts; i++) {
		    auto d = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, i));
            bary += d;
        }

		bary /= nb_verts;

		auto n = normalize(cross(b - a, c - b));
		double dist = dot(p0 - a, n);
		
		if (std::abs(dot(normalize(p0 - a), normalize(bary - a))) < 1e-4)
			continue;

		if (dist < min_dist) {
			min_dist = dist;
			f_idx = mesh_.cells.facet(c_idx, lf);
			lf_idx = lf;
		}

        if (dist < 1e-4)
            break;
	}

	return std::tuple<index_t, index_t>(f_idx, lf_idx);
}

// std::tuple<index_t, index_t> SimpleMeshApplicationExt::pickup_cell_facet(GEO::vec3 p0, index_t c_idx) {
// 	// Search if point is on facet
// 	double min_dist = std::numeric_limits<double>().max();
// 	index_t f_idx = NO_FACET;
// 	index_t lf_idx = NO_FACET;

// 	for (index_t lf = 0; lf < mesh_.cells.nb_facets(c_idx); lf++) {

// 		auto a = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 0));
// 		auto b = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 1));
// 		auto c = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 2));
// 		auto d = mesh_.vertices.point(mesh_.cells.facet_vertex(c_idx, lf, 3));

// 		auto bary = (a + b + c + d) / 4.;

// 		auto n = normalize(cross(b - a, c - b));
// 		double dist = dot(p0 - a, n);
		
// 		if (std::abs(dot(normalize(p0 - a), normalize(bary - a))) < 1e-4)
// 			continue;

// 		if (dist < min_dist) {
// 			min_dist = dist;
// 			f_idx = mesh_.cells.facet(c_idx, lf);
// 			lf_idx = lf;
// 		}
// 	}

// 	return std::tuple<index_t, index_t>(f_idx, lf_idx);
// }

void SimpleMeshApplicationExt::init_rgba_colormap(const std::string& name, int width, int height, unsigned char * data) {
    // like SimpleApplication::init_colormap() but without XPM data, just an u8 array of RGBA values
    colormaps_.push_back(ColormapInfo());
    colormaps_.rbegin()->name = name;
    glGenTextures(1, &colormaps_.rbegin()->texture); // OpenGL will create a new texture and write its id in texture
    glBindTexture(GL_TEXTURE_2D, colormaps_.rbegin()->texture);
    glTexImage2D(
        GL_TEXTURE_2D, 0,
        GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data
    );
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
    );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void SimpleMeshApplicationExt::update_GL_texture(unsigned int texture_index, int width, int height, unsigned char * data) {
    // similar to init_rgba_colormap()
    glBindTexture(GL_TEXTURE_2D, texture_index);
    glTexImage2D(
        GL_TEXTURE_2D, 0,
        GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data
    );
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(
        GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR
    );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
}