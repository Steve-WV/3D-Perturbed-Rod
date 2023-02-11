#include <fstream>
#include <iomanip> 
#include <map>
#include "boost/tokenizer.hpp"

int main(int argc, char *argv[]) {
	
	std::string fname(argv[1]);
	
	std::ifstream is;
	is.open(fname+".node");
	std::string line;
	
	std::ofstream os;
	os.open(fname+".objv");
	os << std::setprecision(18);
	os << "#3d tet data" << std::endl;
	
	std::map<int,int> vmap;
	
	boost::char_separator<char> sep(" ");
	std::getline(is, line); // discard first line
	int v_idx = 0;
	while( std::getline(is, line) ){
		
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line... shouldn't happen, though.
		if ( *ti == "#" ) continue;
		
		int node_idx = std::stoi(*ti);
		vmap[node_idx] = v_idx;

		double x = std::stod(*(++ti));
		double y = std::stod(*(++ti));
		double z = std::stod(*(++ti));
		
		os << "v " << x << " " << y << " " << z << std::endl; 
		
		++v_idx;
	}
	
	is.close();
	
	is.open(fname+".ele");
	std::getline(is, line); // discard first line
	while( std::getline(is, line) ){
		
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line... shouldn't happen, though.
		if ( *ti == "#" ) continue;
		
		int a = vmap[std::stoi(*(++ti))]+1;
		int b = vmap[std::stoi(*(++ti))]+1;
		int c = vmap[std::stoi(*(++ti))]+1;
		int d = vmap[std::stoi(*(++ti))]+1;
		
		os << "t " << a << " " << b << " " << c << " " << d << std::endl; 
		
	}
	
	is.close();
	
	os.close();
	
	
	
}