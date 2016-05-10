#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/IO/print_wavefront.h>
#include "DegradeAnObject.hpp"
#include "TypeDefs.hpp"

int main(int argc, char** argv) {
	srand (time(NULL));
	if(argc != 3) {
		std::cout << "Erreur : ./launch fileToImport nameToExport" << std::endl;
		return 1;
	}
	
	char const *input = argv[1];
	char const *output = argv[2];
	DegradeAnObject o(input, output);
	
	std::vector<Polyhedron> polys = o.getPolyhedrons();
	Polyhedron P = polys[0];
	
	
	Facet f;
	std::cout << o.getFacetFromPoint(P, 1.0, 1.0, 0.5, f) << std::endl;
	std::cout << f.halfedge()->vertex()->point() << std::endl;
	std::cout << f.halfedge()->next()->vertex()->point() << std::endl;
	std::cout << f.halfedge()->next()->next()->vertex()->point() << std::endl;
	std::cout << f.halfedge()->next()->next()->next()->vertex()->point() << std::endl;
	
	
	for(Facet_iterator fi = P.facets_begin() ; fi != P.facets_end() ; ++fi) {
		Point_3 p1 = fi->halfedge()->vertex()->point();
		Point_3 p2 = fi->halfedge()->next()->vertex()->point();
		Point_3 p3 = fi->halfedge()->next()->next()->vertex()->point();
		Point_3 p4 = fi->halfedge()->next()->next()->next()->vertex()->point();
		if(p1 == f.halfedge()->vertex()->point()
		&& p2 == f.halfedge()->next()->vertex()->point() 
		&& p3 == f.halfedge()->next()->next()->vertex()->point() 
		&& p4 == f.halfedge()->next()->next()->next()->vertex()->point()) {
			Halfedge_handle h = P.create_center_vertex(fi->halfedge());
			h->vertex()->point() = Point_3(1.0, 1.0, 0.5);
		}
	}
	
	//h->vertex()->point() = Point_3(h->vertex()->point().x(),h->vertex()->point().y()-1,h->vertex()->point().z());
	
	
	o.exportObj(P);
	
	return 0;
}
