// obj loader from http://jamesgregson.blogspot.fr/2012/05/example-code-for-building.html

#include "TypeDefs.hpp"
#include "DegradeAnObject.hpp"

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

// Constructeur : Charge l'obj et génère les Polyhedron_3 associés
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

// Obj reader, transfère dans les structures de création des polyhedrons les données du fichier
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

// Split un std::string
std::vector<std::string> DegradeAnObject::split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

// Récupère les polyhedrons gréés
std::vector<Polyhedron>& DegradeAnObject::getPolyhedrons() {
	return polys;
}

// Récupère les noms des fichiers
std::vector<std::string> DegradeAnObject::getNames() {
	return names;
}

// Warning : only for triangle mesh. NO. QUAD. MESH.
// Récupère les faces correspondantes au point d'impact et retourne leurs nombres.
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

// Récupère la face de travail du point d'impact
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

// Réalise le travail de raffinage sur une FACE (getFacetsFromPoint == 1)
void DegradeAnObject::refineFacetMesh(Point_3 p, Facet &fs, double epsilon, int index) {
	Facet chkF;
	double maxDistance = distanceBetweenPointAndFacet(p, fs);
	if(maxDistance > epsilon) {
		Halfedge_handle h = splitFacet(fs, index);
		if(getFacetFromPoint(p, chkF, index)) {
			maxDistance = distanceBetweenPointAndFacet(p, chkF);
			if(maxDistance > epsilon) {
				refineFacetMesh(p, chkF, epsilon, index);
			}
			/*else {
				impactAFace(chkF, index, epsilon);
			}
		}
		else {
			impactAFace(fs, index, epsilon);
		}
	}
	else {
		impactAFace(fs, index, epsilon);
	}*/
		}
	}
}

// Subdivise la face (Quaternary subdivision)
Halfedge_handle DegradeAnObject::splitFacet(Facet &fs, int index) {
	Halfedge_handle hh1 = fs.halfedge();
	Halfedge_handle hh2 = fs.halfedge()->next();
	Halfedge_handle hh3 = fs.halfedge()->next()->next();
	std::vector<Halfedge_handle> hhs = splitEdgesOfFacet(fs, index);
	Halfedge_handle h = polys[index].split_facet(hhs[0], hhs[1]);
	h = polys[index].split_facet(h, hhs[2]);
	h = polys[index].split_facet(h, hhs[0]);
	noTVertice(hh1, hh2, hh3, index);
	return h;
}

// Subdivise la face (Barycentric subdivision)
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

// Annule les T Vertices sur les faces adjacentes
void DegradeAnObject::noTVertice(Halfedge_handle hh1, Halfedge_handle hh2, Halfedge_handle hh3, int index) {
	Halfedge_handle h = polys[index].split_facet(hh1->opposite(), hh1->opposite()->next()->next());
	h = polys[index].split_facet(hh2->opposite(), hh2->opposite()->next()->next());
	h = polys[index].split_facet(hh3->opposite(), hh3->opposite()->next()->next());
}

// Casse les arêtes d'une face triangulaire en 2 et les place au centre de son arête. Retourne la liste des pointeurs de ce point
std::vector<Halfedge_handle> DegradeAnObject::splitEdgesOfFacet(Facet fs, int index) {
	Halfedge_handle hh = fs.halfedge();
	std::vector<Halfedge_handle> hhs;
	Point_3 p1 = hh->vertex()->point();
	Point_3 p2 = hh->next()->vertex()->point();
	Point_3 p3 = hh->next()->next()->vertex()->point();
	Halfedge_handle hh1 = polys[index].split_edge(hh);
	hh1->vertex()->point() = meanPoints(p1, p3);
	hhs.push_back(hh1);
	hh = hh->next();
	Halfedge_handle hh2 = polys[index].split_edge(hh);
	hh2->vertex()->point() = meanPoints(p2, p1);
	hhs.push_back(hh2);
	hh = hh->next();
	Halfedge_handle hh3 = polys[index].split_edge(hh);
	hh3->vertex()->point() = meanPoints(p2, p3);
	hhs.push_back(hh3);
	return hhs;
}

// Calcule une moyenne de deux points
Point_3 DegradeAnObject::meanPoints(Point_3 p1, Point_3 p2) {
	Point_3 pt(p2.x() + p1.x(), p2.y() + p1.y(), p2.z() + p1.z());
	pt = Point_3(pt.x()/2, pt.y()/2, pt.z()/2);
	return pt;
}

// Calcule une moyenne d'un vecteur de plusieurs points
Point_3 DegradeAnObject::meanPoints(std::vector<Point_3> points) {
	Point_3 pt(0.0, 0.0, 0.0);
	int size = points.size();
	for(int i = 0 ; i < size ; i++) {
		pt = Point_3(pt.x() + points[i].x(), pt.y() + points[i].y(), pt.z() + points[i].z());
	}
	pt = Point_3(pt.x()/size, pt.y()/size, pt.z()/size);
	return pt;
}

// Evalue la distance entre un point et les sommets d'une face. Retourne la distance avec le sommet le plus éloigné
double DegradeAnObject::distanceBetweenPointAndFacet(Point_3 p, Facet f) {
	double maxDistance = 0;
	double distance;
	std::vector<Point_3> pts = getAllPointsFromFacet(f);
	for(int i = 0 ; i < pts.size() ; i++) {
		distance = CGAL::squared_distance(p, pts[i]);
		if(maxDistance < distance) {
			maxDistance = distance;
		}
	}
	return sqrt(maxDistance);
}

// Réalise l'impact d'une face
void DegradeAnObject::impactAFace(Facet &fs, int index, double epsilon) {
	Vector_3 normal = getNormalOfFacet(fs);
	Halfedge_handle h = barycentricMesh(fs, index);
	Point_3 pt = h->vertex()->point();
	h->vertex()->point() = Point_3(pt.x() - 10*normal.x(), pt.y() - 10*normal.y(), pt.z() - 10*normal.z());
}

// Retourne le vecteur normal d'une face
Vector_3 DegradeAnObject::getNormalOfFacet(Facet fs) {
	Vector_3 normal = CGAL::cross_product(fs.halfedge()->next()->vertex()->point() - fs.halfedge()->vertex()->point(), fs.halfedge()->next()->next()->vertex()->point() - fs.halfedge()->next()->vertex()->point());
	return normal;
}

// Normalise un vecteur
Vector_3 DegradeAnObject::normalizeVector(Vector_3 v) {
	double norm = sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
	return v/norm;
}

// Génère des points autour du point d'impact prévu
std::vector<Point_3> DegradeAnObject::generatePointsOnFacet(Point_3 p, double ray, Facet fs) {
	std::vector<Point_3> pts;
	pts.push_back(p);
	Vector_3 normal = normalizeVector(getNormalOfFacet(fs));
	Kernel::Plane_3 pl(fs.halfedge()->vertex()->point(), normal);
	Vector_3 orth = (pl.base1()) * ray;
	Point_3 chkPt = p + orth;
	int i = 0;
	Facet test;
	double teta = 10;
	double c = cos(teta);
	double s = sin(teta);
	Kernel::RT m00 = normal.x() * normal.x() * (1-c) + c;
	Kernel::RT m01 = normal.x() * normal.y() * (1-c) - normal.z() * s;
	Kernel::RT m02 = normal.x() * normal.z() * (1-c) + normal.y() * s;
	Kernel::RT m10 = normal.x() * normal.y() * (1-c) + normal.z() * s;
	Kernel::RT m11 = normal.y() * normal.y() * (1-c) + c;
	Kernel::RT m12 = normal.z() * normal.y() * (1-c) - normal.x() * s;
	Kernel::RT m20 = normal.x() * normal.z() * (1-c) - normal.y() * s;
	Kernel::RT m21 = normal.z() * normal.y() * (1-c) + normal.x() * s;
	Kernel::RT m22 = normal.z() * normal.z() * (1-c) + c;
	CGAL::Aff_transformation_3<Kernel> rotate(m00, m01, m02, m10, m11, m12, m20, m21, m22);
	while(i < 36) {
		orth = rotate.transform(orth);
		chkPt = p + orth;
		//std::cout << chkPt << std::endl;
		pts.push_back(chkPt);
		i++;
	}
	return pts;
}

void DegradeAnObject::impactTheFacetArea(Point_3 p, Facet fs, double ray, int index) {
	double distance = 0;
	Vector_3 normal = normalizeVector(getNormalOfFacet(fs));
	Kernel::Plane_3 pl(fs.halfedge()->vertex()->point(), normal);
	for(Point_iterator pi = polys[index].points_begin() ; pi != polys[index].points_end() ; ++pi) {
		Point_3 chkPt = *pi;
		if(pl.has_on(chkPt)) {
			distance = sqrt(CGAL::squared_distance(p, chkPt));
			if(distance < ray) {
				*pi = Point_3(chkPt.x() - 0.1*normal.x(), chkPt.y() - 0.1*normal.y(), chkPt.z() - 0.1*normal.z());
			}
		}
	}
}

// Récupère la liste de tous les points d'une face
std::vector<Point_3> DegradeAnObject::getAllPointsFromFacet(Facet f) {
	std::vector<Point_3> pts;
	Halfedge_handle hh = f.halfedge();
	pts.push_back(hh->vertex()->point());
	hh = hh->next();
	while(hh->vertex()->point() != pts[0]) {
		pts.push_back(hh->vertex()->point());
		hh = hh->next();
	}
	return pts;
}

double DegradeAnObject::rangeRandomAlg2 (int min, int max) {
    int n = max - min + 1;
    int remainder = RAND_MAX % n;
    int x;
    do {
        x = rand();
    }
    while(x >= RAND_MAX - remainder);
    return min + x % n;
}

// Start the process
void DegradeAnObject::startDeformation() {
	std::vector<Facet> fs;
	std::vector<int> indexes;
	Point_3 p(1.0, 0.562, 0.818);
	double ray = 0.1;
	getFacetsFromPoint(p, fs, indexes);
	
	if(fs.size() == 1) { // On a facet
		std::vector<Point_3> pts = generatePointsOnFacet(p, ray, fs[0]);
		for(int i = 0 ; i < pts.size() ; i++) {
			refineFacetMesh(pts[i], fs[0], 0.02, indexes[0]);
		}
		impactTheFacetArea(p, fs[0], ray, indexes[0]);
	}
}

// Exporte les polyhedrons en obj. REVOIR UN DETAIL DESSUS VIS A VIS DE LA NUMEROTATION
void DegradeAnObject::exportObj() {
	std::ofstream ofs(output);
	for(int i = 0 ; i < polys.size() ; i++) {
		CGAL::print_polyhedron_wavefront(ofs, polys[i]);
	}
	ofs.close();
	std::cout << "objet exporté" << std::endl;
	
}

// Exporte un polyhedron ciblé en obj. Utile à des fins de tests
void DegradeAnObject::exportObj(Polyhedron P) {
	std::ofstream ofs(output);
	CGAL::print_polyhedron_wavefront(ofs, P);
	ofs.close();
	std::cout << "objet exporté" << std::endl;
	
}
