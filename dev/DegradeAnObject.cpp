#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <math.h>
 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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

void DegradeAnObject::exportObj(Polyhedron P) {
	std::ofstream ofs(output);
	CGAL::print_polyhedron_wavefront(ofs, P);
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
	if(myFile.is_open()) {
		readFromThisObject(myFile, "", 0);
		myFile.close();
	}
	std::reverse(minFacets.begin(), minFacets.end());
	std::reverse(names.begin(), names.end());
	std::reverse(coords.begin(), coords.end());
	std::reverse(faces.begin(), faces.end());
}

void DegradeAnObject::readFromThisObject(std::ifstream & myFile, std::string n, int it) {
	int curMinFacet = -1;
	std::string line;
	std::vector<double> vertexes;
	std::vector< std::vector<int> > facets;
	std::string name = n;

	while(getline(myFile,line)){
		if(line.size() >= 2) {
			if(line.at(0) == 'v' && (line.at(1) != 't' && line.at(1) != 'n')){
				std::vector<std::string> vert = split(line, ' ');
				vertexes.push_back(::atof(vert[1].c_str()));
				vertexes.push_back(::atof(vert[2].c_str()));
				vertexes.push_back(::atof(vert[3].c_str()));
			}
			else if(line.at(0) == 'f'){
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
			else if(line.at(0) == 'o' || line.at(0) == 'g') {
				readFromThisObject(myFile, line.substr(2), it+1);
			}	
		}
	}
	if(it != 0) {
		minFacets.push_back(curMinFacet);
		names.push_back(name);
		coords.push_back(vertexes);
		faces.push_back(facets);
	}
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
/* obsolete
char* DegradeAnObject::subchars(char* str, short x, short y){
    char* ret = new char[y+1];
    for(short i=x; i<x+y; i++)
        ret[i-x]=str[i];
    ret[y] = '\0';
    return ret;
}
*/
std::vector<Polyhedron>& DegradeAnObject::getPolyhedrons() {
	return polys;
}

std::vector<std::string> DegradeAnObject::getNames() {
	return names;
}

// Warning : only for triangle mesh. NO. QUAD. MESH.
int DegradeAnObject::getFacetsFromPoint(Point_3 p, std::vector<Facet> &fs, std::vector<int> &index) {
	for(int i = 0 ; i < polys.size() ; i++) {
		for(Facet_iterator fi = polys[i].facets_begin(); fi != polys[i].facets_end() ; ++fi) {
			Point_3 p1 = fi->halfedge()->vertex()->point();
			Point_3 p2 = fi->halfedge()->next()->vertex()->point();
			Point_3 p3 = fi->halfedge()->next()->next()->vertex()->point();
			Kernel::Triangle_3 t3(p1, p2, p3);
			if(t3.has_on(p)) {
				index.push_back(i);
				fs.push_back(*fi);
			}
		}
	}
	return fs.size();
}

int DegradeAnObject::getFacetFromPoint(Point_3 p, Facet &fs, int index) {
	for(Facet_iterator fi = polys[index].facets_begin(); fi != polys[index].facets_end() ; ++fi) {
		Point_3 p1 = fi->halfedge()->vertex()->point();
		Point_3 p2 = fi->halfedge()->next()->vertex()->point();
		Point_3 p3 = fi->halfedge()->next()->next()->vertex()->point();
		Kernel::Triangle_3 t3(p1, p2, p3);
		if(t3.has_on(p)) {
			fs = *fi;
			return 1;
		}
	}
	return 0;
}

void DegradeAnObject::refineFacetMesh(Point_3 p, Facet &fs, double epsilon, int index) {
	Facet chkF;
	Halfedge_handle h = splitFacet(fs, index);
	std::cout << getFacetFromPoint(p, chkF, index) << std::endl;
	if(distanceBetweenPointAndFacet(p, chkF.halfedge()->vertex()->point()) > epsilon) {
		refineFacetMesh(p, chkF, epsilon, index);
	}
	else {
		std::cout << "stop" << std::endl;
		//impactAFace(p, chkF, index);
	}
}

Halfedge_handle DegradeAnObject::splitFacet(Facet fs, int index) {
	splitEdgesOfFacet(fs, index);
	Halfedge_handle h = barycentricMesh(fs, index);
	//noTVertice(fs, index);
	return h;
}

Halfedge_handle DegradeAnObject::barycentricMesh(Facet fs, int index) {
	std::vector<Point_3> points;
	Halfedge_handle hh = fs.halfedge();
	Point_3 p1 = hh->vertex()->point();
	points.push_back(p1);
	hh = hh->next();
	while(hh->vertex()->point() != p1) {
		points.push_back(hh->vertex()->point());
		hh = hh->next();
	}
	Halfedge_handle h = polys[index].create_center_vertex(fs.halfedge());
	h->vertex()->point() = meanPoints(points);
	return h;
}

void DegradeAnObject::noTVertice(Facet fs, int index) {
	std::cout << fs.halfedge()->vertex()->point() << std::endl;
	barycentricMesh(*(fs.halfedge()->opposite()->facet()), index);
	//barycentricMesh(*(fs.halfedge()->next()->opposite()->facet()), index);
	//barycentricMesh(*(fs.halfedge()->next()->next()->opposite()->facet()), index);
}

void DegradeAnObject::splitEdgesOfFacet(Facet fs, int index) {
	Halfedge_handle hh = fs.halfedge();
	Point_3 p1 = hh->vertex()->point();
	Point_3 p2 = hh->next()->vertex()->point();
	Point_3 p3 = hh->next()->next()->vertex()->point();
	Halfedge_handle hh1 = polys[index].split_edge(hh);
	hh1->vertex()->point() = meanPoints(p1, p3);
	hh = hh->next();
	Halfedge_handle hh2 = polys[index].split_edge(hh);
	hh2->vertex()->point() = meanPoints(p2, p1);
	hh = hh->next();
	Halfedge_handle hh3 = polys[index].split_edge(hh);
	hh3->vertex()->point() = meanPoints(p2, p3);
}

Point_3 DegradeAnObject::meanPoints(Point_3 p1, Point_3 p2) {
	Point_3 pt(p2.x() + p1.x(), p2.y() + p1.y(), p2.z() + p1.z());
	pt = Point_3(pt.x()/2, pt.y()/2, pt.z()/2);
	return pt;
}

Point_3 DegradeAnObject::meanPoints(std::vector<Point_3> points) {
	Point_3 pt(0.0, 0.0, 0.0);
	int size = points.size();
	for(int i = 0 ; i < size ; i++) {
		pt = Point_3(pt.x() + points[i].x(), pt.y() + points[i].y(), pt.z() + points[i].z());
	}
	pt = Point_3(pt.x()/size, pt.y()/size, pt.z()/size);
	return pt;
}

double DegradeAnObject::distanceBetweenPointAndFacet(Point_3 p, Point_3 pfs) {
	return sqrt((p.x() - pfs.x()) * (p.x() - pfs.x()) + (p.y() - pfs.y()) * (p.y() - pfs.y()) + (p.z() - pfs.z()) * (p.z() - pfs.z()));
}

void DegradeAnObject::impactAFace(Point_3 p, Facet &fs, int index) {
	Halfedge_handle h = polys[index].create_center_vertex(fs.halfedge());
	h->vertex()->point() = Point_3	(1.0-0.1, 0.5, 0.8);;
}

void DegradeAnObject::changeAllPoints() {
	for(std::vector<Polyhedron>::iterator P = polys.begin() ; P != polys.end() ; ++P) {
		for ( Vertex_iterator v = P->vertices_begin(); v != P->vertices_end(); ++v) {
			Point_3 p(v->point().x()+((double) rand() / (RAND_MAX)),v->point().y()+((double) rand() / (RAND_MAX)),v->point().z()+((double) rand() / (RAND_MAX)));
			v->point() = p;
			std::cout << v->point() << std::endl;
		}
		std::cout << "--" << std::endl;
	}
}

