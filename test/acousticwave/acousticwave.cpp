#include <models/acoustics/acoustics.hpp>
#include <models/acoustics/fluxes/godunov.hpp>
#include <solvers/finitevolume.hpp>
#include <mesh/io/meshreader.hpp>

#include <iostream>
#include <string>
#include <ctime>
#include <cmath>

using namespace std;
using namespace ConservationLaw2D;

typedef double real_t;
typedef Model::LinearAcoustics<real_t>			myModel;
typedef NumericalFlux::Godunov<myModel>			myNumFlux;
typedef Solver::FiniteVolume<myModel,myNumFlux>	mySolver;
typedef mySolver::FVMesh						myMesh;

// Tipo per la soluzione (vettore d-dim)
typedef myModel::SolType	SolType;

// Stato iniziale in variabili (p,u,v)
inline SolType init( size_t color, real_t x, real_t y ) {
	SolType sol = SolType::Zero();
	// p, u, v
	sol[0] = 2.0*exp(-80.0*(x*x+y*y));
	sol[1] = 0.0;
	sol[2] = 0.0;
	return sol;
}

// Condizioni al bordo
inline SolType bc( SolType& wl, size_t color, real_t x, real_t y, real_t nx, real_t ny, real_t t ) {
	return wl;
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
	// Definisco il modello con parametri RHO0 e K0
	myModel model(1.0,1.0);
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
	// Inizializzo il solutore
	solver.init();
	solver.setDirectory("./data");
	// Passi temporali
	for (int i = 0; i < 1000; ++i) {
		std::cout << "== Timestep " << i << " == currtime: " << std::setw(8) << solver.getCurrTime();
		std::cout << ", dt = " << std::setw(8) << solver.getCurrDt() << std::endl;
		solver.timestep();
		if (i%50 == 0) solver.framegrab(i/50, gnuplot, interpolated);
	}
	return EXIT_SUCCESS;
}
