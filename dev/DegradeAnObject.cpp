// obj loader from http://jamesgregson.blogspot.fr/2012/05/example-code-for-building.html

#include "TypeDefs.hpp"
#include "DegradeAnObject.hpp"


/**
 * Methodes pour charger l'obj dans les structures de données
 * 
 */
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

/**
 * Méthodes de raffinage de maillage
 * 
 */

// Warning : only for triangle mesh. NO. QUAD. MESH.
// Récupère les faces correspondantes au point d'impact et retourne leurs nombres.
int DegradeAnObject::getFacetsFromPoint(Point_3 p, std::vector<Facet> &fs, std::vector<int> &index) {
	for(int i = 0 ; i < polys.size() ; i++) {
		for(Facet_iterator fi = polys[i].facets_begin(); fi != polys[i].facets_end() ; ++fi) {
			Point_3 p1 = fi->halfedge()->vertex()->point();
			Point_3 p2 = fi->halfedge()->next()->vertex()->point();
			Point_3 p3 = fi->halfedge()->next()->next()->vertex()->point();
			if(isAPointInATriangle(p, p1, p2, p3)) {
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
		if(isAPointInATriangle(p, p1, p2, p3)) { 
			fs = *fi;
			return 1;
		}
	}
	return 0;
}

// Place un point p sur la face fs, et relie p aux sommets de fs.
Halfedge_handle DegradeAnObject::putAPointOnAFacet(Point_3 p, int index) {
	Facet fs;
	getFacetFromPoint(p, fs, index);
	Halfedge_handle h = polys[index].create_center_vertex(fs.halfedge());
	h->vertex()->point() = p;
	return h;
}

// Place un point p sur la face fs, et relie p aux sommets de fs.
bool DegradeAnObject::isAPointOnThisFacet(Point_3 p, Facet fs, int index) {
	std::vector<Point_3> pts = getAllPointsFromFacet(fs);
	if(isAPointInATriangle(p, pts[0], pts[1], pts[2])) {
		return true;
	}
	else {
		return false;
	}
}

bool DegradeAnObject::sameSide(Point_3 p1, Point_3 p2, Point_3 a, Point_3 b) {
    Vector_3 cp1 = CGAL::cross_product(b-a, p1-a);
    Vector_3 cp2 = CGAL::cross_product(b-a, p2-a);
    if(cp1*cp2 >= 0) {
		return true;
	}
    else {
		return false;
	}
}

bool DegradeAnObject::isAPointInATriangle(Point_3 p, Point_3 a, Point_3 b, Point_3 c) {
	Plane_3 pl(a, b, c);
	if(pl.has_on(p)) {
		if(sameSide(p, a, b, c) && sameSide(p, b, a, c) && sameSide(p, c, a, b)) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}
}

bool DegradeAnObject::twoPointsOnTheFacet(Point_3 p1, Point_3 p2, Facet fs, int index) {
	return (isAPointOnThisFacet(p1, fs, index) && isAPointOnThisFacet(p2, fs, index));
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

Halfedge_handle DegradeAnObject::splitEdge(Halfedge_handle hh, Point_3 p, int index) {
	Halfedge_handle rh = polys[index].split_edge(hh);
	rh->vertex()->point() = p;
	return rh;
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
		distance = squared_distance(p, pts[i]);
		if(maxDistance < distance) {
			maxDistance = distance;
		}
	}
	return sqrt(maxDistance);
}

// Distance entre 2 points
double DegradeAnObject::squared_distance(Point_3 p1, Point_3 p2) {
	return to_double((p1.x() - p2.x()) * (p1.x() - p2.x()) + (p1.y() - p2.y()) * (p1.y() - p2.y()) + (p1.z() - p2.z()) * (p1.z() - p2.z()));
}

// Retourne le vecteur normal d'une face
Vector_3 DegradeAnObject::getNormalOfFacet(Facet fs) {
	Vector_3 normal = CGAL::cross_product(fs.halfedge()->next()->vertex()->point() - fs.halfedge()->vertex()->point(), fs.halfedge()->next()->next()->vertex()->point() - fs.halfedge()->next()->vertex()->point());
	return normal;
}

// Normalise un vecteur
Vector_3 DegradeAnObject::normalizeVector(Vector_3 v) {
	double norm = sqrt(to_double(v.x() * v.x() + v.y() * v.y() + v.z() * v.z()));
	return v/norm;
}

// Génère des points autour du point d'impact prévu
std::vector<Point_3> DegradeAnObject::generatePointsOnFacet(Point_3 p, double ray, Facet fs, int nbPts) {
	std::vector<Point_3> pts;
	Vector_3 normal = normalizeVector(getNormalOfFacet(fs));
	Kernel::Plane_3 pl(fs.halfedge()->vertex()->point(), normal);
	Vector_3 orth = (pl.base1()) * ray;
	Vector_3 smallOrth;
	Point_3 chkPt;
	chkPt = p + orth;
	pts.push_back(chkPt);
	Facet test;
	double teta = M_PI/(nbPts/2.0);
	for(int i = 1 ; i < nbPts ; i++) {
		orth = rotationVector(orth, normal, teta);
		chkPt = p + orth;
		pts.push_back(chkPt);
	}
	return pts;
}

// Réalise un impact sur la face fs à partir d'une liste de points à déplacer, répartis par couronne (pts[0] = première couronne intérieure, pts[0][0] = premier point de la première couronne)
void DegradeAnObject::impactTheFacetArea(std::vector< std::vector<Point_3> > pts, Facet fs, double ray, int index) {
	double str = 0.02;
	Vector_3 normal = normalizeVector(getNormalOfFacet(fs));
	Kernel::Plane_3 pl(fs.halfedge()->vertex()->point(), normal);
	for(int i = 0 ; i < pts.size() ; i++) {
		for(int j = 0 ; j < pts[i].size() ; j++) {
			bool chk = false;
			Point_iterator pi = polys[index].points_begin();
			while(!chk) {
				++pi;
				if(*pi == pts[i][j]) {
					*pi = Point_3(pi->x() - (impactStrengh(str, i))*normal.x(), pi->y() - (impactStrengh(str, i))*normal.y(), pi->z() - (impactStrengh(str, i))*normal.z());
					chk = true;
				}
			}
		}
	}
}

double DegradeAnObject::impactStrengh(double initStrengh, int i) {
	return initStrengh+initStrengh*(i/(i+1.0));
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

// Génère un nombre aléatoire entre min et max
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

// Génère le vecteur qui a subi une rotation d'angle teta
Vector_3 DegradeAnObject::rotationVector(Vector_3 v, Vector_3 normal, double teta) {
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
	return rotate.transform(v);
}

// Dessine tous les points qui déterminent l'impact à réaliser, en plusieurs couronnes.
void DegradeAnObject::drawImpactOnFacet(Point_3 p, double ray, std::vector<Point_3> pts, Facet initFs, int index, int nbCouronnes) {
	std::vector< std::vector<Point_3> > points;
	std::vector<Point_3> tmp;
	std::vector<Halfedge_handle> hhs;
	Halfedge_handle hh;
	hh = putAPointOnAFacet(pts[0], index);
	tmp.push_back(hh->vertex()->point());
	for(int i = 1 ; i < pts.size() ; i++) {
		hh = joinTwoPoints(pts[i], pts[i-1], index, tmp);
	}
	joinFirstAndLast(pts[0], pts[pts.size()-1], index, tmp);
	for(int i = 1 ; i < nbCouronnes ; i++) {
		hhs = getHalfedgesOfPoints(tmp, index);
		tmp = drawInsideImpactOnFacet(tmp, hhs, initFs, index);
		points.push_back(tmp);
	}
	impactTheFacetArea(points, initFs, ray, index);
}

// Dessine la couronne intérieure
std::vector<Point_3> DegradeAnObject::drawInsideImpactOnFacet(std::vector<Point_3> points, std::vector<Halfedge_handle> hhs, Facet f, int index) {
	std::vector<Point_3> pts;
	for(int i = 0 ; i < points.size() ; i++) {
		int j;
		if(i == points.size()-1) {
			j = 0;
		}
		else {
			j = i+1;
		}
		Vector_3 h(hhs[i]->opposite()->vertex()->point(), hhs[i]->vertex()->point());
		Vector_3 g(hhs[j]->opposite()->vertex()->point(), hhs[j]->vertex()->point());
		Vector_3 norm = getNormalOfFacet(f);
		Vector_3 rh = normalizeVector(rotationVector(h, norm, M_PI/2));
		Vector_3 rg = normalizeVector(rotationVector(g, norm, M_PI/2));
		Vector_3 comb = 0.01*normalizeVector(rh+rg);
		Point_3 newPoint = hhs[i]->vertex()->point() + comb;
		Halfedge_handle hh = polys[index].split_vertex(hhs[j]->opposite(), hhs[i]);
		hh->vertex()->point() = newPoint;
		polys[index].split_facet(hh->opposite()->next()->next(), hh->opposite());
		polys[index].split_facet(hh->next()->next(), hh);
		pts.push_back(newPoint);
	}
	return pts;
}

// Relie 2 points dans un polyèdre. p2 est le point PRECEDENT à p1.
Halfedge_handle DegradeAnObject::joinTwoPoints(Point_3 p1, Point_3 p2, int index, std::vector<Point_3> & pts) {
	Halfedge_handle hh;
	std::vector<Facet> fcts;
	std::vector<int> indexes;
	Halfedge_handle prevHalf;
	Facet fs;
	getFacetFromPoint(p1, fs, index);
	if(twoPointsOnTheFacet(p1, p2, fs, index)) {
		prevHalf = putAPointOnAFacet(p1, index);
		pts.push_back(prevHalf->vertex()->point());
	}
	else {
		fcts.clear();
		Segment_3 s(p1, p2);
		getFacetsFromPoint(p2, fcts, indexes);
		hh = getExteriorHalfedge(p2, s, fcts);
		Halfedge_handle previousHalfedge = hh->next();
		Halfedge_handle newEdge = addAndJoinNewPoint(p2, previousHalfedge, hh, s, index);
		pts.push_back(newEdge->vertex()->point());
		prevHalf = joinTwoPoints(p1, newEdge->vertex()->point(), index, pts);
	}
	
	return prevHalf;
}

// Relie le premier et le dernier point
Halfedge_handle DegradeAnObject::joinFirstAndLast(Point_3 p1, Point_3 p2, int index, std::vector<Point_3> & pts) {
	Halfedge_handle hh;
	bool chk = false;
	std::vector<Facet> fcts;
	std::vector<int> indexes;
	Halfedge_handle prevHalf;
	Facet fs;
	getFacetsFromPoint(p1, fcts, indexes);
	for(int i = 0 ; i < fcts.size() ; i++) {
		if(twoPointsOnTheFacet(p1, p2, fcts[i], index)) {
			chk = true;
		}
	}
	if(!chk) {
		fcts.clear();
		Segment_3 s(p1, p2);
		getFacetsFromPoint(p2, fcts, indexes);
		hh = getExteriorHalfedge(p2, s, fcts);
		Halfedge_handle previousHalfedge = hh->next();
		Halfedge_handle newEdge = addAndJoinNewPoint(p2, previousHalfedge, hh, s, index);
		pts.push_back(newEdge->vertex()->point());
		prevHalf = joinFirstAndLast(p1, newEdge->vertex()->point(), index, pts);
	}
	
	return prevHalf;
}

// Recherche les halfedges des - facets du point - qui ne contiennent pas le point
Halfedge_handle DegradeAnObject::getExteriorHalfedge(Point_3 p, Segment_3 s, std::vector<Facet> fcts) {
	Halfedge_handle retHh;
	for(int i = 0 ; i < fcts.size() ; i++) {
		Halfedge_handle hh = fcts[i].halfedge();
		for(int j = 0 ; j < 3 ; j++) {
			if(hh->vertex()->point() != p && hh->opposite()->vertex()->point() != p) {
				Segment_3 seg(hh->opposite()->vertex()->point(), hh->vertex()->point());
				if(!seg.is_degenerate()) {
					if(CGAL::do_intersect(s, seg)) {
						retHh = hh;
					}
				}
			}
			hh = hh->next();
		}
	}
	return retHh;
}

Halfedge_handle DegradeAnObject::addAndJoinNewPoint(Point_3 p, Halfedge_handle previousHalfedge, Halfedge_handle hh, Segment_3 s, int index) {
	Point_3 intersect;
	Halfedge_handle splittedHalfedge;
	Segment_3 seg(hh->opposite()->vertex()->point(), hh->vertex()->point());
	Point_3* chkPt; 
	CGAL::cpp11::result_of<Kernel::Intersect_3(Segment_3, Segment_3)>::type result = CGAL::intersection(s, seg);
	if (result) {
		chkPt = boost::get<Point_3 >(&*result);
		intersect = *chkPt;
	}
	Halfedge_handle split = splitEdge(hh, intersect, index);
	Halfedge_handle hhx = polys[index].split_facet(previousHalfedge, split);
	Halfedge_handle oppositePoint = hhx->next()->opposite();
	polys[index].split_facet(oppositePoint, oppositePoint->next()->next());
	
	return oppositePoint;
}

std::vector<Halfedge_handle> DegradeAnObject::getHalfedgesOfPoints(std::vector<Point_3> points, int index) {
	std::vector<Halfedge_handle> hhs;
	Halfedge_handle hh = getHalfedgeBetweenTwoPoints(points[0], points[points.size()-1], index);
	hhs.push_back(hh);
	for(int i = 1 ; i < points.size() ; i++) {
		hhs.push_back(getHalfedgeBetweenTwoPoints(points[i], points[i-1], index));
	}
	return hhs;
}

Halfedge_handle DegradeAnObject::getHalfedgeBetweenTwoPoints(Point_3 p1, Point_3 p2, int index) {
	bool found = false;
	Halfedge_iterator hi = polys[index].halfedges_begin();
	while(!found) {
		if(p1 == hi->vertex()->point() && p2 == hi->opposite()->vertex()->point()) {
			found = true;
		}
		else {
			hi++;
		}
	}
	return hi;
}

/**
 * Public methods
 * 
 */

// Starts the process
void DegradeAnObject::startDeformation() {
	std::vector<Facet> fs;
	std::vector<int> indexes;
	Point_3 p(1.0, 0.562, 0.818);
	double ray = 0.1;
	int nbCouronnes = 4;
	int nbPoints = 10;
	getFacetsFromPoint(p, fs, indexes);
	if(fs.size() == 1) { // On a facet
		Facet f = fs[0];
		std::vector<Point_3> pts = generatePointsOnFacet(p, ray, f, nbPoints);
		drawImpactOnFacet(p, ray, pts, f, indexes[0], nbCouronnes);
	}
}

// Exporte les polyhedrons en obj. REVOIR UN DETAIL DESSUS VIS A VIS DE LA NUMEROTATION
void DegradeAnObject::exportObj() {
	std::ofstream ofs(output);
	for(int i = 0 ; i < polys.size() ; i++) {
		ofs << "o " << names[i] << "\n";
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
