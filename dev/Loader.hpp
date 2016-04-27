#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/IO/print_wavefront.h>

class ObjToPolyhedron {
	
	typedef CGAL::Simple_cartesian<double>     Kernel;
	typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
	typedef Polyhedron::HalfedgeDS             HalfedgeDS;
	
	Polyhedron P;
	std::vector<double> coords;
	std::vector< std::vector<int> >    faces;
	char const *output;
	
	
	public : 
	ObjToPolyhedron(char const *input, char const *output);	
	int get_first_integer( const char *v );
	void load_obj( const char *filename, std::vector<double> &coords, std::vector< std::vector<int> > &faces );
	void exportObj();
	
};
