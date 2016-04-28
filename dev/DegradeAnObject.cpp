#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
 
#include<CGAL/Simple_cartesian.h>
#include<CGAL/Polyhedron_incremental_builder_3.h>
#include<CGAL/Polyhedron_3.h>
#include<CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/IO/print_wavefront.h>
#include "Loader.hpp"

#include "TypeDefs.hpp"
 
// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:
	std::vector<double> &coords;
	std::vector< std::vector<int> >    &faces;
	polyhedron_builder( std::vector<double> &_coords, std::vector< std::vector<int> > &_faces ) : coords(_coords), faces(_faces) {}
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
				B.add_vertex_to_facet( faces[i][j] );
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
	load_obj( input, coords, faces );
	if( coords.size() == 0 ) {
		std::cout << "Aucun objet n'a été chargé" << std::endl;
	}
	else {
		std::cout << "objet chargé" << std::endl;
	  
	 // build a polyhedron from the loaded arrays

		polyhedron_builder<HalfedgeDS> builder( coords, faces );
		P.delegate( builder );
	}
}

void DegradeAnObject::exportObj() {
	std::ofstream ofs(output);
	CGAL::print_polyhedron_wavefront(ofs, P);
	ofs.close();
	std::cout << "objet exporté" << std::endl;
	
	/* EXPORT OFF
	// write the polyhedron out as a .OFF file
	std::ofstream os("dump.off");
	os << P;
	os.close();*/
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
void DegradeAnObject::load_obj( const char *filename, std::vector<double> &coords, std::vector< std::vector<int> > &faces ){
	double x, y, z;
	char line[1024], str[1024];

	// open the file, return if open fails
	FILE *fp = fopen(filename, "r" );
	if( !fp ) return;
  
	// read lines from the file, if the first character of the
	// line is 'v', we are reading a vertex, otherwise, if the
	// first character is a 'f' we are reading a facet
	while( fgets( line, 1024, fp ) ){
		if(line[0] == "o") {
			char *token = std::strtok(line, " ");
			while (token != NULL) {
				if(token[0] != 'o') {
					names.push_back(token);
				}
				token = std::strtok(NULL, " ");
			}
		}
		else if(line[0] == 'v' && (line[1] != 't' && line[1] != 'n')){
			sscanf( line, "%*s%lf%lf%lf", &x, &y, &z );
			coords.push_back( x );
			coords.push_back( y );
			coords.push_back( z );
		}
		else if( line[0] == 'f' ){
			char *token = std::strtok(line, " ");
			std::vector<int> f;
			while (token != NULL) {
				if(token[0] != 'f') {
					f.push_back(get_first_integer(token)-1);
				}
				token = std::strtok(NULL, " ");
			}
			faces.push_back(f);
			/*
			sscanf( line, "%*s%s%s%s%s", v0, v1, v2, v3 );
			faces.push_back( get_first_integer( v0 )-1 );
			faces.push_back( get_first_integer( v1 )-1 );
			faces.push_back( get_first_integer( v2 )-1 );
			faces.push_back( get_first_integer( v3 )-1 );
			*/
		}
	}
	fclose(fp); 
}

Polyhedron DegradeAnObject::getPolyhedron() {
	return P;
}

void DegradeAnObject::changeAllPoints() {
	for ( Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v) {
		Point_3 p(v->point().x()+((double) rand() / (RAND_MAX)),v->point().y()+((double) rand() / (RAND_MAX)),v->point().z()+((double) rand() / (RAND_MAX)));
        v->point() = p;
        std::cout << v->point() << std::endl;
	}
}
