#include "TypeDefs.hpp"

class DegradeAnObject {
	
	std::vector<std::string> names;
	std::vector<Polyhedron> polys;
	std::vector< std::vector<double> > coords;
	std::vector< std::vector< std::vector<int> > > faces; // for each polyhedron, for each face, each coord.
	char const *output;
	std::vector<int> minFacets;
	
	void load_obj( const char *filename );
	void readFromThisObject(std::ifstream & myFile, std::string n, int it);
	std::vector<std::string> split(const std::string &s, char delim);
	
	int getFacetsFromPoint(Point_3 p, std::vector<Facet> &fs, std::vector<int> &index);
	int getFacetFromPoint(Point_3 p, Facet &fs, int index);
	Halfedge_handle putAPointOnAFacet(Point_3 p, int index);
	bool isAPointOnThisFacet(Point_3 p, Facet fs, int index);
	bool sameSide(Point_3 p1, Point_3 p2, Point_3 a, Point_3 b);
	bool isAPointInATriangle(Point_3 p, Point_3 a, Point_3 b, Point_3 c);
	bool twoPointsOnTheFacet(Point_3 p1, Point_3 p2, Facet fs, int index);
	void noTVertice(Halfedge_handle hh1, Halfedge_handle hh2, Halfedge_handle hh3, int index);
	std::vector<Halfedge_handle> splitEdgesOfFacet(Facet fs, int index);
	Halfedge_handle splitEdge(Halfedge_handle hh, Point_3 p, int index);
	Point_3 meanPoints(Point_3 p1, Point_3 p2);
	Point_3 meanPoints(std::vector<Point_3> points);
	double distanceBetweenPointAndFacet(Point_3 p, Facet f);
	double squared_distance(Point_3 p1, Point_3 p2);
	Vector_3 getNormalOfFacet(Facet fs);
	Vector_3 normalizeVector(Vector_3 v);
	std::vector<Point_3> generatePointsOnFacet(Point_3 p, double ray, Facet fs, int nbPts);
	void impactTheFacetArea(std::vector< std::vector<Point_3> > pts, Facet fs, double ray, int index);
	double impactStrengh(double initStrengh, int i);
	std::vector<Point_3> getAllPointsFromFacet(Facet f);
	double rangeRandomAlg2 (int min, int max);
	Vector_3 rotationVector(Vector_3 v, Vector_3 normal, double teta);
	void drawImpactOnFacet(Point_3 p, double ray, std::vector<Point_3> pts, Facet initFs, int index, int nbCouronnes); 
	std::vector<Point_3> drawInsideImpactOnFacet(std::vector<Point_3> points, std::vector<Halfedge_handle> hhs, Facet f, int index);
	Halfedge_handle joinTwoPoints(Point_3 p1, Point_3 p2, int index, std::vector<Point_3> & pts);
	Halfedge_handle joinFirstAndLast(Point_3 p1, Point_3 p2, int index, std::vector<Point_3> & pts);
	Halfedge_handle getExteriorHalfedge(Point_3 p, Segment_3 s, std::vector<Facet> fcts);
	Halfedge_handle addAndJoinNewPoint(Point_3 p, Halfedge_handle previousHalfedge, Halfedge_handle hh, Segment_3 s, int index);
	std::vector<Halfedge_handle> getHalfedgesOfPoints(std::vector<Point_3> points, int index);
	Halfedge_handle getHalfedgeBetweenTwoPoints(Point_3 p1, Point_3 p2, int index);
	
	public : 
	DegradeAnObject(char const *input, char const *output);	
	void exportObj();
	void exportObj(Polyhedron P);
	void startDeformation();
};
