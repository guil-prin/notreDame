#include "TypeDefs.hpp"
#include "DegradeAnObject.hpp"

int main(int argc, char** argv) {
	srand (time(NULL));
	if(argc != 3) {
		std::cout << "Erreur : ./launch fileToImport nameToExport" << std::endl;
		return 1;
	}
	
	char const *input = argv[1];
	char const *output = argv[2];
	DegradeAnObject o(input, output);

	o.startDeformation();
	/*else if(fs.size() == 2) {
		for(int i = 0 ; i < polys.size() ; i++) {
			Halfedge_handle hh;
			for(Halfedge_iterator hi = polys[i].halfedge_begin() ; hi != polys[i].halfedge_begin_end() ; ++fi) {
			}
			Halfedge_handle hh = polys[i].split_edge(fi->halfedge());
			hh->vertex()->point() = Point_3(1.0-0.2, 1.0-0.2, 0.5);
		}
	}*/
	
	//h->vertex()->point() = Point_3(h->vertex()->point().x(),h->vertex()->point().y()-1,h->vertex()->point().z());
	
	
	o.exportObj();
	
	return 0;
}
