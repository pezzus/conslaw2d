#include <models/shallowwater/shallowwater.hpp>
#include <solvers/fluxes/laxfriedrichs.hpp>
#include <solvers/finitevolume.hpp>
#include <mesh/io/meshreader.hpp>

#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;
using namespace ConservationLaw2D;

typedef double real_t;
typedef Model::ShallowWater<real_t>              myModel;
typedef NumericalFlux::LaxFriedrichs<myModel>    myNumFlux;
typedef Solver::FiniteVolume<myModel,myNumFlux>  mySolver;
typedef mySolver::FVMesh                         myMesh;

// Tipo per la soluzione (vettore d-dim)
typedef myModel::SolType	SolType;

// Stato iniziale in variabili (rho,u,v)
inline SolType init( size_t color, real_t x, real_t y ) {
	SolType sol = SolType::Zero();
	// rho, u, v
	sol[0] = ( x < 0 ) ? 4.0 : 1.0;
	sol[1] = 0.0;
	sol[2] = 0.0;
	return sol;
}

// Condizioni al bordo
inline SolType bc( SolType& wl, size_t color, real_t x, real_t y, real_t nx, real_t ny, real_t t ) {
	SolType wr = SolType::Zero();
	switch (color) {
		case 8:
		case 11:
		case 15:
		case 10:
		case 18:
		case 13:
			// Riflettente
			wr[0] = wl[0];
			wr[1] = (ny*ny-nx*nx)*wl[1] - 2*nx*ny*wl[2];
			wr[2] = - 2*nx*ny*wl[1] + (nx*nx-ny*ny)*wl[2];
			break;
		default:
			// Assorbente
			wr = wl;
			break;
	}
	return wr;
}

int main(int argc, char **argv) {
	// Parametri in ingresso
	bool gnuplot(false), interpolated(false);
	string meshfile;
	if ( argc < 2 ) {
		cout << "Usage: " << argv[0] << " [options] meshfile.msh" << endl;
		cout << "Options:" << endl;
		cout << "  --gnuplot\t\tGenerate plot and animation from gnuplot" << endl;
		cout << "  --interpolated\t\tInterpolate solution on vertices" << endl;
		exit(1);
	}
	for (int i=1; i<argc; ++i) {
		if (!strcmp(argv[i],"--gnuplot")) 
			gnuplot = true;
		else if (!strcmp(argv[i],"--interpolated"))
			interpolated = true;
		else
			meshfile = argv[i];
	}
	// Definisco il modello
	myModel model;
	// Definisco la mesh
	myMesh mesh;
	// Leggo la mesh
	Mesh::IO::MeshReader(mesh, meshfile);
	// Definisco il solutore per il mio modello
	mySolver solver(model, mesh);
	// Inizializzo alcuni parametri
	solver.setCFLmax(0.1);
	solver.setIC(init);
	solver.setBC(bc);
	solver.init();
	solver.setDirectory("./data");
	// Passi temporali
	for (int i = 0; i < 1000; ++i) {
		std::cout << "== Timestep " << i << " == currtime: " << std::setw(8) << solver.getCurrTime();
		std::cout << ", dt = " << std::setw(8) << solver.getCurrDt() << std::endl;
		solver.timestep();
		if (i%50 == 0) solver.framegrab(i/50, gnuplot, interpolated);
	}
	return 0;
}
