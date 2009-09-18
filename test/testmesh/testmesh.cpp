#include <mesh/mesh_default_traits.hpp>
#include <mesh/io/meshreader.hpp>

#include <iostream>
#include <ctime>

using namespace std;
using namespace ConservationLaw2D;

// Caratterizzazione degli Elementi della Mesh
typedef Mesh::DefaultTraits<double>	Traits;
typedef Traits::PolygonalMesh		SimpleMesh;

int main(int argc, char **argv) {
	// Avvertimento
	if ( argc < 2 ) {
		cout << "Usage: " << argv[0] << " meshfile.msh" << endl;
		exit(EXIT_SUCCESS);
	}
	// Cronometro
	clock_t ck0, ck1;
	// Leggo la mesh
	cout << "====== READING MESH ======" << endl;
	SimpleMesh m;
	ck0=clock();
	Mesh::IO::MeshReader(m, argv[1]);
	ck1=clock();
	m.stats();
	cout << "Time: " << double(ck1-ck0)/CLOCKS_PER_SEC << " seconds" << endl;
	return 0;
}
