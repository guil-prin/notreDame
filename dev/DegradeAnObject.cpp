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

#include "DegradeAnObject.hpp"
#include "TypeDefs.hpp"
 
// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
	std::vector<double> &coords;
	std::vector< std::vector<int> >    &faces;
	int minFacet;
	polyhedron_builder( std::vector<double> &_coords, std::vector< std::vector<int> > &_faces, int _minFacet ) : coords(_coords), faces(_faces), minFacet(_minFacet) {}
	void operator()( HDS& hds) {
		typedef typename HDS::Vertex   Vertex;
		typedef typename Vertex::Point Point;

		// create a cgal incremental builder
		CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
		B.begin_surface( coords.size()/3, faces.size() );

		// add the polyhedron vertices
		for( int i=0; i<(int)coords.size(); i+=3 ){
			B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
		}
		// add the polyhedron faces
		for(int i = 0 ; i < faces.size() ; i++) {
			B.begin_facet();
			for(int j = 0 ; j < faces[i].size() ; j++) {
				B.add_vertex_to_facet( faces[i][j] - minFacet );
			}
			B.end_facet();
		}

		// finish up the surface
		B.end_surface();
	}
};

DegradeAnObject::DegradeAnObject(char const *input, char const *output) {
	
	this->output = output;
	// load the input file
	load_obj(input);
	if( coords.size() == 0 ) {
		std::cout << "Aucun objet n'a été chargé" << std::endl;
	}
	else {
		std::cout << "objet chargé" << std::endl;
		std::reverse(minFacets.begin(), minFacets.end());
		std::reverse(names.begin(), names.end());
		std::reverse(coords.begin(), coords.end());
		std::reverse(faces.begin(), faces.end());
	 // build polyhedrons from the loaded arrays
		for(int i = 0 ; i < names.size() ; i++) {
			Polyhedron P;
			polyhedron_builder<HalfedgeDS> builder( coords[i], faces[i], minFacets[i] );
			P.delegate( builder );
			polys.push_back(P);
		}
	}
}

void DegradeAnObject::exportObj() {
	std::ofstream ofs(output);
	for(int i = 0 ; i < polys.size() ; i++) {
		CGAL::print_polyhedron_wavefront(ofs, polys[i]);
	}
	ofs.close();
	std::cout << "objet exporté" << std::endl;
	
}
 
// reads the first integer from a string in the form
// "334/455/234" by stripping forward slashes and
// scanning the result
int DegradeAnObject::get_first_integer( const char *v ){
	 int ival;
	 std::string s( v );
	 std::replace( s.begin(), s.end(), '/', ' ' );
	 sscanf( s.c_str(), "%d", &ival );
	 return ival;
}
 
// barebones .OFF file reader, throws away texture coordinates, normals, etc.
// stores results in input coords array, packed [x0,y0,z0,x1,y1,z1,...] and
// faces array packed [T0a,T0b,T0c,T1a,T1b,T1c,...]
void DegradeAnObject::load_obj(const char *filename){
	std::string line;
	std::ifstream myFile (filename);
	if(myFile.is_open())
	{
		readFromThisObject(myFile);
		myFile.close();
	}
}

void DegradeAnObject::readFromThisObject(std::ifstream & myFile) {
	int curMinFacet = -1;
	std::string line;
	std::vector<double> vertexes;
	std::vector< std::vector<int> > facets;
	std::string name;

	while(getline(myFile,line)){
		if(line.at(0) == 'v' && (line.at(1) != 't' && line.at(1) != 'n')){
			std::vector<std::string> vert = split(line, ' ');
			vertexes.push_back(::atof(vert[1].c_str()));
			vertexes.push_back(::atof(vert[2].c_str()));
			vertexes.push_back(::atof(vert[3].c_str()));
		}
		else if( line.at(0) == 'f' ){
			std::vector<std::string> face = split(line, ' ');
			std::vector<int> f;
			for(int i = 1 ; i < face.size() ; i++) {
				int j = ::atoi(face[i].c_str()) - 1;
				f.push_back(j);
				if(curMinFacet == -1) {
					curMinFacet = j;
				}
				else if(j < curMinFacet) {
					curMinFacet = j;
				}
			}
			facets.push_back(f);
		}
		else if(line.at(0) == 'o') {
			name = line.substr(2);
			readFromThisObject(myFile);
		}	
	}
	minFacets.push_back(curMinFacet);
	names.push_back(name);
	coords.push_back(vertexes);
	faces.push_back(facets);
}

std::vector<std::string> DegradeAnObject::split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

char* DegradeAnObject::subchars(char* str, short x, short y){
    char* ret = new char[y+1];
    for(short i=x; i<x+y; i++)
        ret[i-x]=str[i];
    ret[y] = '\0';
    return ret;
}

/*
Polyhedron DegradeAnObject::getPolyhedrons() {
	return polys;
}
*/
std::vector<std::string> DegradeAnObject::getNames() {
	return names;
}

/*
void DegradeAnObject::changeAllPoints() {
	for ( Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v) {
		Point_3 p(v->point().x()+((double) rand() / (RAND_MAX)),v->point().y()+((double) rand() / (RAND_MAX)),v->point().z()+((double) rand() / (RAND_MAX)));
        v->point() = p;
        std::cout << v->point() << std::endl;
	}
}
*/
