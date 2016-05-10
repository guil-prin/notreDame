#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/print_wavefront.h>

#include "TypeDefs.hpp"

class DegradeAnObject {
	
	std::vector<std::string> names;
	std::vector<Polyhedron> polys;
	std::vector< std::vector<double> > coords;
	std::vector< std::vector< std::vector<int> > > faces; // for each polyhedron, for each face, each coord.
	char const *output;
	std::vector<int> minFacets;
	
	
	public : 
	DegradeAnObject(char const *input, char const *output);	
	int get_first_integer( const char *v );
	void load_obj( const char *filename );
	void readFromThisObject(std::ifstream & myFile, std::string n, int it);
	std::vector<std::string> split(const std::string &s, char delim);
	void exportObj();
	void exportObj(Polyhedron P);
	
	std::vector<Polyhedron> getPolyhedrons();
	std::vector<std::string> getNames();
	
	bool getFacetFromPoint(Polyhedron P, double x, double y, double z, Facet &f);
	void changeAllPoints();
	//char* subchars(char* str, short x, short y);
};
