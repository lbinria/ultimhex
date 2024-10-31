#pragma once

#include "mesh_metadata.h"

struct Tool {

	/**
	 * Check whether the tool is compatible with current mesh according to its metadata
	 */
	virtual bool is_compatible(MeshMetadata mesh_metadata);

	/**
	 * Draw GUI
	 */
	virtual void draw_gui(/* pass mesh, tet, hex ... */);

	/**
	 * Call when user press escape key
	 */
	virtual void escape_callback();

	virtual void mouse_button_callback(int button, int action, int mods, int source);
	virtual void hover_callback(double x, double y, int source);
	virtual void key_callback(int key, int scancode, int action, int mods);

};