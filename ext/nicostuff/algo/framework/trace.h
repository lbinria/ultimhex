#ifndef __TRACE_H__
#define __TRACE_H__


#include "basic.h"



#include <ultimaille/all.h>
#include <iostream>
#include <string>

#include <fstream>
#include <math.h>

#if __cplusplus >= 199711L  
#include <filesystem>
#endif

using namespace UM;


namespace Trace {
	extern bool drop_mesh_is_active ;
	extern bool trace_steps_active ;
	extern int num_drop ;
	extern std::string outputdir;


	inline std::string lua_filename() {return outputdir + "/view.lua";}

	inline void abort(std::string msg= "Abort called" ) { throw(std::runtime_error(msg.c_str())); }
	//    _      _ _   _      _ _         _   _          
	//   (_)_ _ (_) |_(_)__ _| (_)_____ _| |_(_)___ _ __ 
	//   | | ' \| |  _| / _` | | |_ / _` |  _| / _ \ '  | 
	//   |_|_||_|_|\__|_\__,_|_|_/__\__,_|\__|_\___/_||_|
	void initialize(std::string p_outputdir, std::string p_outputprefix = "");

	//    _                _               _   ___ ___ 
	//   | |   ___  __ _  | |_ _____ __   /_\ | _ \_ _|
	//   | |__/ _ \/ _` | |  _/ -_) \ /  / _ \|  _/| | 
	//   |____\___/\__, |  \__\___/_\_\ /_/ \_\_| |___|
	//             |___/                              
	struct SwitchTextInScope {
		SwitchTextInScope(bool make_active);
		~SwitchTextInScope();
		static void force_activity(bool b); // dirty, doesn't stack change
		bool save_val;
	};


	enum ForceDropActive {
		DROP_UNCHANGED,
		DROP_ON,
		DROP_OFF,
	};
	struct Section {
		Section(const std::string& str, ForceDropActive drop = DROP_UNCHANGED);
		~Section();
		bool drop_save_val;
		bool trace_was_active;
	};
	void step(std::string stepname, int alert_level = 0);
	void alert(std::string msg);

	void log_value(std::string const& str, double val, int alert_level = 0);
	void log_string(std::string const& str, std::string const& val, int alert_level = 0);
	void show_log();
	void append_py_log(std::string const& filename);
	void append_log(std::string const& filename);


	//    ___                 __  __        _               _   ___ ___ 
	//   |   \ _ _ ___ _ __  |  \/  |___ __| |_  ___ ___   /_\ | _ \_ _|
	//   | |) | '_/ _ \ '_ \ | |\/| / -_|_-< ' \/ -_|_-<  / _ \|  _/| | 
	//   |___/|_| \___/ .__/ |_|  |_\___/__/_||_\___/__/ /_/ \_\_| |___|
	//                |_|                                               
	struct SwitchDropInScope {
		SwitchDropInScope(bool make_active);
		~SwitchDropInScope();
		bool save_val;
		static void force_activity(bool b); // dirty, doesn't stack change
	};






	


};


#endif