#include <iostream> // std::cout, std::cin 
#include <ostream>  // << 
#include <istream> 
#include <fstream>
#include<sstream>
#include <fstream>
#include <iomanip> 
#include <map>
#include "boost/tokenizer.hpp"
using namespace std; 

int main(int argc, char *argv[]) {
	
	std::string fname(argv[1]);
	
	std::ifstream is;
	is.open(fname+".mesh");
	//is.open(fname+".msh");
	//is.open("test.mesh");
	std::string line;
	
	

	std::ofstream os;
	os.open(fname+".objv");
	//os.open("test.objv");	
	os << std::setprecision(18);
	os << "#3d tet data" << std::endl;
	
	std::map<int,int> vmap;
	
	boost::char_separator<char> sep(" ");
	std::getline(is, line); // discard four lines
	std::getline(is, line); 
	std::getline(is, line); 
	std::getline(is, line);
	std::getline(is, line);
	int v_idx = 0;
	while( std::getline(is, line) ){
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();

				
		if (ti==te) continue; //empty line... shouldn't happen, though.
		if ( *ti == "Triangles" ) break; // we're done with vertices...
		
		//int node_idx = std::stoi(*ti);
		//vmap[node_idx] = v_idx;
	
		double x = std::stod(*(ti++));
		double y = std::stod(*(ti++));
	std::cout << "Here" << std::endl;
		double z = std::stod(*(ti++));
		
		
		os << "v " << x << " " << y << " " << z << std::endl; 

		
		
		++v_idx;
	}

	std::cout << "Finished" << std::endl;
	
	while( std::getline(is, line) ){
		
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line... shouldn't happen, though.
		if ( *ti == "Tetrahedra" ) break; // we're done with triangles...
	}
	
	std::getline(is, line); // discard one line
	while( std::getline(is, line) ){
		
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line... shouldn't happen, though.
		if ( *ti == "End" ) break;
		
		int a = std::stoi(*(ti++));
		int b = std::stoi(*(ti++));
		int c = std::stoi(*(ti++));
		int d = std::stoi(*(ti++));
		
		os << "t " << a << " " << b << " " << c << " " << d << std::endl; 
		
	}
	
	is.close();
	
	os.close();
	
	
	
}
