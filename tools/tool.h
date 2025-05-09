#pragma once

#include "../context.h"
#include "../mesh_metadata.h"

struct Tool {

	virtual std::string get_name() = 0;

	Tool(Context &ctx) : ctx(ctx) {}

	/**
	 * Check whether the tool is compatible with current mesh according to its metadata
	 */
	virtual bool is_compatible() = 0;

	/**
	 * Draw GUI
	 * Return true if selected tool has changed
	 */
	virtual bool draw_object_properties() = 0;

	/**
	 * Draw viewer property for this tool
	 */
	virtual void draw_viewer_properties() = 0;

	/**
	* Draw GL
	*/
	virtual void draw(GEO::vec4f hovered_color, GEO::vec4f selected_color, GEO::SimpleApplication::ColormapInfo colorMapInfo) = 0;

	/**
	 * Call when user press escape key
	 */
	virtual void escape_callback() = 0;

	virtual void validate_callback() = 0;

	virtual void mouse_button_callback(int button, int action, int mods, int source) = 0;
	virtual void scroll_callback(double xoffset, double yoffset) = 0;
	virtual void hover_callback(double x, double y, int source) = 0;
	virtual void key_callback(int key, int scancode, int action, int mods) = 0;

	/**
	 * Clear tool
	 */
	virtual void clear() = 0;

	Context &ctx;

};