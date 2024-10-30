#include <framework/trace.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <chrono>
#include <ctime> 
using namespace UM;



namespace Trace {

	bool trace_steps_active = true;
	int num_drop = 0;
	std::string outputdir = "";
	std::string outputprefix = "";
	bool drop_mesh_is_active = false;

	struct  LogTime{
		/**
		 * CheckPoint is the list item of LogTime
		 * "up" is it's father
		 * "right" is the next item at the same level
		 */
		struct CheckPoint {
			CheckPoint(std::string const& p_n, unsigned int p_up) { n = p_n; up = p_up; t = clock(); right = (unsigned int)(-1); }
			std::string n;
			clock_t t;
			unsigned int right;
			unsigned int up;
		};

		static std::vector<CheckPoint> check;
		static std::vector<std::pair<std::string, double> > out_values;
		static std::vector<std::pair<std::string, std::string> > out_strings;


		static bool is_start_section(unsigned int i){
		return check[i].right != i + 1;
	}

		static bool is_end_section(unsigned int i){
		return check[i].n == "end section";
	}
		static bool is_final(unsigned int i){
		return i + 1 == check.size();
	}

		static double time(unsigned int i){
		return (double(check[check[i].right].t) - double(check[i].t)) / double(CLOCKS_PER_SEC);
	}

		static unsigned int dec(unsigned int i = (unsigned int)(-1)) {
		if (i == (unsigned int)(-1))
			i = (unsigned int)(check.size() - 1);
		unsigned int res = 0;
		i = check[i].up;
		while (i != (unsigned int)(-1)) {
			res++;
			i = check[i].up;
		}
		return res;
	}

		static unsigned int lastdec() {
		if (!check.empty())
			return dec((unsigned int)(check.size() - 1));
		return 0;
	}

		static void debug() {
		std::cerr << std::endl;
		std::cerr << "--------------BEGIN DEBUG-------------------" << std::endl;
		for (size_t i = 0; i < check.size(); i++) {
			std::cerr << std::string((unsigned int)(4 * dec((unsigned int)(i))), ' ') << i << "  r = " << check[i].right << " u = " << check[i].up
				<< "\tstart" << check[i].t << "\tname" << check[i].n << std::endl;
		}
		std::cerr << std::endl;
		std::cerr << "-------------- END  DEBUG-------------------" << std::endl;
	}

		static std::string cur_stack() {
		std::string res;
		std::vector<unsigned int> stack;
		{
			unsigned int i = (unsigned int)(check.size() - 1);
			while (i != (unsigned int)(-1)) { stack.push_back(i); i = check[i].up; }
		}
		for (int i = int(stack.size()) - 1; i >= 0; i--) {
			res.append(check[stack[size_t(i)]].n);
			if (i > 0) res.append(" ==> ");
		}
		return res;
	}

	static void report(std::ostream& out, unsigned int timing_depth=10000) {
		if (check.empty()) return;
		if (check.back().n != "the end") { add_step("the end"); }

		if (timing_depth != (unsigned int)(-1)) {
			out << "\n***********************************************************" << std::endl;
			out << "                  TIMING SUMMARY " << std::endl;
			for (unsigned int i = 0; i < check.size() - 1; i++) {
				if (dec(i) > timing_depth) continue;
				if (is_start_section(i))
					out << std::string(4U * dec(i), ' ') << time(i) << "\t====>  " << check[i].n << std::endl;
				else if (!is_end_section(i) && check[i].n != "begin section")
					out << std::string(4U * dec(i), ' ') << time(i) << "\t" << check[i].n << std::endl;
			}				
			double ns = double(check[check.size() - 1U].t - check[0].t) / double(CLOCKS_PER_SEC);
			log_value("Time", ns);
			out << ns << "\tTOTAL" << std::endl;
		}
		out << "\n***********************************************************" << std::endl;
		out << "                  OUPUT VALUES" << std::endl;
		for (auto const& v : out_values)
			out << v.second << " \t" << v.first << std::endl;
	}



	static void report_py(std::ostream& out, unsigned int timing_depth) {
		if (!check.empty() && check.back().n != "the end") { add_step("the end"); }
		out << "{\"finished\": \"yes\"";
		if (timing_depth != (unsigned int)(-1)) {
			for (unsigned int i = 0; i < check.size() - 1; i++) {
				if (dec(i) > timing_depth) continue;
				if (is_start_section(i)) out << ",\"TIME_" << check[i].n << "\" :  " << time(i);
				else if (!is_end_section(i) && check[i].n != "begin section")
					out << ",\"TIME_" << check[i].n << "\":  " << time(i);
			}
		}
		for (auto const& v : out_values)
			out << ", \"" << v.first << "\": " << v.second;
		for (auto const& v : out_strings)
			out << ",\"" << v.first << "\": \"" << v.second << "\"";
		out << "}";

	}

	static void report_log(std::ostream& out, unsigned int timing_depth) {
		if (check.empty()) {out << "finished=yes\n"; return; }
		if (!check.empty() && check.back().n != "the end") { add_step("the end"); }
		out << "finished=yes\n";
		if (timing_depth != (unsigned int)(-1)) {
			for (unsigned int i = 0; i < check.size() - 1; i++) {
				if (dec(i) > timing_depth) continue;
				if (is_start_section(i)) out << "TIME_" << check[i].n << "=" << time(i)<<std::endl;
				else if (!is_end_section(i) && check[i].n != "begin section")
					out << "TIME_" << check[i].n << "=" << time(i) << std::endl;
			}
		}
		for (auto const& v : out_values)
			out << v.first << "=" << v.second << std::endl;
		for (auto const& v : out_strings)
			out << v.first << "=" << v.second <<std::endl;
		

	}





	// construct API
		static void log_value(std::string const& str, double val) {
		std::cerr << "LogTime >> " << str << " = " << val << std::endl;
		out_values.push_back(std::pair<std::string, double>(str, val));
	}
		static void log_string(std::string const& str, std::string const& val) {
		std::cerr << "LogTime >> " << str << " = " << val << std::endl;
		for (auto& record : out_strings) if (record.first.compare(str) == 0) {
			record.second = record.second + "," + val;
			return;
		}
		out_strings.push_back(std::pair<std::string, std::string>(str, val));
	}

		static void add_step(std::string const& name) {
		if (check.empty()) {
			CheckPoint c(name, (unsigned int)(-1));
			check.push_back(c);
		} else {
			CheckPoint c(name, check.back().up);
			check.back().right = (unsigned int)(check.size());
			check.push_back(c);
		}
		const char* symbol = "#>=_-~..............";
		std::cerr << std::endl << std::string(4 * lastdec(), ' ') << std::string(80 - 4 * lastdec(), symbol[dec()]) << std::endl;
		std::cerr << std::string(4 * lastdec() + 4, ' ') << cur_stack() << std::endl;
		auto glo = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		std::cerr << std::string(4 * lastdec() + 4, ' ') << std::ctime(&glo) << std::endl;
	}

		static void start_section(std::string const& secname, std::string const& name="begin section") {
		add_step(secname);
		CheckPoint c(name, (unsigned int)(check.size()) - 1);
		check.push_back(c);
	}
		static void end_section() {
		unsigned int u = check.back().up;
		check[u].right = (unsigned int)(check.size());
		check.back().right = (unsigned int)(check.size());
		CheckPoint c("end section", check[u].up);
		check.push_back(c);
	}

	
		static void drop_file(std::string const& filename, bool append, unsigned int timing_depth) {
		std::ofstream f;
		if (append) f.open(filename.c_str(), std::fstream::app);
		else		f.open(filename.c_str());
		report(f, timing_depth);
		f.close();
	}
	};

	std::vector<LogTime::CheckPoint> LogTime::check;
	std::vector<std::pair<std::string, double> > LogTime::out_values;
	std::vector<std::pair<std::string, std::string> > LogTime::out_strings;








	void SwitchDropInScope::force_activity(bool b) { drop_mesh_is_active = b; }

	SwitchDropInScope::SwitchDropInScope(bool make_active) {
		if (outputdir == "") make_active = false;
		save_val = drop_mesh_is_active;
		drop_mesh_is_active = make_active;
	};
	SwitchDropInScope::~SwitchDropInScope() {
		drop_mesh_is_active = save_val;
	}

	void SwitchTextInScope::force_activity(bool b) { trace_steps_active = b; }
	SwitchTextInScope::SwitchTextInScope(bool make_active) {
		if (outputdir == "") make_active = false;
		save_val = trace_steps_active;
		trace_steps_active = make_active;
	};
	SwitchTextInScope::~SwitchTextInScope() {
		trace_steps_active = save_val;
	}



	void initialize(std::string p_outputdir, std::string p_outputprefix) {
		outputdir = p_outputdir;
		drop_mesh_is_active = true;
		trace_steps_active = true;
		std::string lua_file_name = outputdir + std::string("/view.lua");
		std::string dos_compatible;
		for (auto c : lua_file_name) if (c == '/') dos_compatible += "\\"; else  dos_compatible += c;
		system(("del " + dos_compatible).c_str());
	};



	void set_alert_level(int alert_level) {
		if (alert_level == 0) shell_color(10, 0);
		if (alert_level == 1) shell_color(14, 0);
		if (alert_level == 2) shell_color(12, 0);
		if (alert_level == 3) shell_color(0, 12);
	}
	void step(std::string stepname, int alert_level) {
		set_alert_level(alert_level);
		LogTime::add_step(stepname);
		shell_color(15, 0);
	}
	Section::Section(const std::string& str,ForceDropActive drop_mode) {
		drop_save_val = drop_mesh_is_active;
		if (drop_mode != DROP_UNCHANGED) {
			drop_mesh_is_active = (drop_mode == DROP_ON);
		}
		set_alert_level(0);
		LogTime::start_section(str);
		shell_color(15, 0);
	}
	Section::~Section() {
		
		drop_mesh_is_active = drop_save_val;
		LogTime::end_section(); 
	}

	void log_value(std::string const& str, double val, int alert_level) { set_alert_level(alert_level); LogTime::log_value(str, val); shell_color(15, 0);}
	void log_string(std::string const& str, std::string const& val, int alert_level) {set_alert_level(alert_level);LogTime::log_string(str, val); shell_color(15, 0);}
	void alert(std::string msg) { log_string("ALERT",msg,2); }
	void show_log() {
		LogTime::report(std::cerr);
	}
	void append_py_log(std::string const& filename) {
		std::ofstream myfile;
		myfile.open(filename, std::ios::out);
		LogTime::report_py(myfile, -1);
		myfile.close();
	}
	void append_log(std::string const& filename) {
		std::ofstream myfile;
		myfile.open(filename, std::ios::out);
		LogTime::report_log(myfile, -1);
		myfile.close();
	}









};