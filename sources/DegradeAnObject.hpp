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
	
	std::vector<Polyhedron>& getPolyhedrons();
	std::vector<std::string> getNames();
	
	int getFacetsFromPoint(Point_3 p, std::vector<Facet> &fs, std::vector<int> &index);
	int getFacetFromPoint(Point_3 p, Facet &fs, int index);
	void refineFacetMesh(Point_3 p, Facet &fs, double epsilon, int index);
	Halfedge_handle splitFacet(Facet fs, int index);
	Halfedge_handle barycentricMesh(Facet fs, int index);
	void noTVertice(Halfedge_handle hh1, Halfedge_handle hh2, Halfedge_handle hh3, int index);
	void splitEdgesOfFacet(Facet fs, int index);
	Point_3 meanPoints(Point_3 p1, Point_3 p2);
	Point_3 meanPoints(std::vector<Point_3> points);
	double distanceBetweenPointAndFacet(Point_3 p, Point_3 pfs);
	void impactAFace(Facet &fs, int index);
	void changeAllPoints();
	//char* subchars(char* str, short x, short y);
};
