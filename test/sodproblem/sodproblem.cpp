#include <models/eulero/eulero.hpp>
#include <solvers/fluxes/laxfriedrichs.hpp>
#include <models/eulero/fluxes/godunovROE.hpp>
#include <models/eulero/fluxes/godunovHLL.hpp>
#include <models/eulero/fluxes/godunovHLLC.hpp>
#include <models/eulero/fluxes/rusanov.hpp>
#include <solvers/finitevolume.hpp>
#include <mesh/io/meshreader.hpp>

#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;
using namespace ConservationLaw2D;

typedef double real_t;
typedef Model::Eulero<real_t>					myModel;
//typedef NumericalFlux::Rusanov<myModel>		myNumFlux;
//typedef NumericalFlux::LaxFriedrichs<myModel>	myNumFlux;
typedef NumericalFlux::GodunovRoe<myModel>		myNumFlux;
//typedef NumericalFlux::GodunovHLL<myModel>	myNumFlux;
//typedef NumericalFlux::GodunovHLLC<myModel>	myNumFlux;
typedef Solver::FiniteVolume<myModel,myNumFlux>	mySolver;
typedef mySolver::FVMesh						myMesh;

// Tipo per la soluzione (vettore d-dim)
typedef myModel::SolType	SolType;

// Stato iniziale in variabili (p,u,v)
inline SolType init( size_t color, real_t x, real_t y ) {
	SolType sol;
	// rho, u, v
	sol[0] = (x < 0) ? 1.0 : 0.125;
	sol[1] = (x < 0) ? 0.75 : 0.0;
	sol[2] = 0.0;
	sol[3] = (x < 0) ? 1.0 : 0.1;
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
	// Definisco il modello (con gamma=1.4)
	myModel model(1.4);
	// Definisco la mesh
	myMesh mesh;
	// Leggo la mesh
	Mesh::IO::MeshReader(mesh, meshfile);
	// Storing della geometria e output di alcune statistiche
	mesh.init_geom();
	mesh.stats();
	// Definisco il solutore per il mio modello
	mySolver solver(model, mesh);
	// Inizializzo alcuni parametri
	solver.setCFLmax(0.1);
	solver.setIC(init);
	solver.setBC(bc);
	solver.setDirectory("./data");
	// Inizializzo il solutore
	solver.init();
	// Passi temporali
	for (int i = 0; i <= 500; ++i) {
		std::cout << "== Timestep " << i << " == currtime: " << std::setw(8) << solver.getCurrTime();
		std::cout << ", dt = " << std::setw(8) << solver.getCurrDt() << std::endl;
		solver.timestep();
		if (i%50 == 0) solver.framegrab(i/50, gnuplot, interpolated);
	}
	return 0;
}

