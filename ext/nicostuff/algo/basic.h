#ifndef __BASIC_H__
#define __BASIC_H__

#include <ultimaille/all.h>

#include <string>
#include<vector>
#include<set>
#include<map>


#ifdef WIN32
#include <windows.h>
inline void shell_color(int t, int f) {
	HANDLE H = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(H, WORD(f * 16 + t));
}
#undef min
#undef max
#else
inline void shell_color(int t, int f) {
	if (t == 15 && f == 0) {
		std::cerr << "\e[0m";
		return;
	}
	std::cerr << "\e[0";
	if (t & 8) std::cerr << ";1";
	t = ((t & 1) << 2) | (t & 2) | ((t >> 2) & 1);
	std::cerr << ";3" << char('0' + char(t));
	f = ((f & 1) << 2) | (f & 2) | ((f >> 2) & 1);
	if (f != 0) std::cerr << ";4" << char('0' + char(f));
	std::cerr << 'm';
}
#endif


inline double rand_range(double min, double max) { return min + (max - min) * double(rand() % 10000) / 10000.; }

#define FOR(i, n) for(int i = 0; i < (int) n; i++)

inline std::ostream& operator<<(std::ostream& os, const UM::LinExpr& line){
	FOR(i,line.size()) {
		os<<" ";
		if (i>0 && line[i].value>0)
		os<<"+";
		os<<line[i].value;
		if (line[i].index!=-1) os<<"X_"<<line[i].index;
		
	} return os;
}

// easy debug output by "plop"
template <class T> inline void plop_value(const T& val) { std::cerr << val; }
template <class T> inline void plop_value(const std::vector<T>& data) { std::cerr << "\n"; FOR(i, data.size()) std::cerr << "[" << i << "] -> " << data[i] << "\n"; }
template <class T> inline void plop_value(const std::set<T>& data) { std::cerr << "\n"; for (auto e : data) std::cerr << e << "\n"; }
template <class T, class Q> inline void plop_value(const std::map<T, Q>& data) { std::cerr << "\n"; for (auto e : data) std::cerr << "key= " << e.first << "  =val=>  " << e.second << "\n"; }

inline std::string plop_unquote_string(const char* str) {
	std::string result(str);
	if (result.length() > 2 && result[0] == '\"' && result[result.length() - 1] == '\"') {
		result = result.substr(1, result.length() - 2);
	}
	return result;
}
#define plop(x) {shell_color(14,0); std::cerr << "|plop|=>"<<" line:"<< __LINE__ <<"  "<< plop_unquote_string(#x) <<" : ";plop_value(x);std::cerr << "     in  "<< __FILE__ <<std::endl;shell_color(15,0);}





using namespace UM;

template  <class T>
struct ToolBox {
	ToolBox(T& a) {}
	static void NoToolboxDefinedForThisObject() { um_assert(!"Please do not call this IDE helper"); }
};





#endif