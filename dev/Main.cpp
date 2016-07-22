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
	
	o.exportObj();
	
	return 0;
}
