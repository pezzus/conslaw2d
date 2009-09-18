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
//typedef NumericalFlux::Rusanov<myModel>			myNumFlux;
//typedef NumericalFlux::LaxFriedrichs<myModel>		myNumFlux;
typedef NumericalFlux::GodunovRoe<myModel>			myNumFlux;
//typedef NumericalFlux::GodunovHLL<myModel>		myNumFlux;
//typedef NumericalFlux::GodunovHLLC<myModel>		myNumFlux;
typedef Solver::FiniteVolume<myModel,myNumFlux>	mySolver;
typedef mySolver::FVMesh						myMesh;

// Tipo per la soluzione (vettore d-dim)
typedef myModel::SolType	SolType;

// Stato iniziale in variabili (p,u,v)
inline SolType init( size_t color, real_t x, real_t y ) {
	SolType sol;
	// rho, u, v
	sol[0] = 1.4;
	sol[1] = 2.9;
	sol[2] = 0.0;
	sol[3] = 1.0;
	return sol;
}

// Condizioni al bordo
inline SolType bc( SolType& wl, size_t color, real_t x, real_t y, real_t nx, real_t ny, real_t t ) {
	SolType wr;
	switch (color) {
		case 0:
			// Riflettente
			wr[0] = wl[0];
			wr[1] = (ny*ny-nx*nx)*wl[1] - 2*nx*ny*wl[2];
			wr[2] = - 2*nx*ny*wl[1] + (nx*nx-ny*ny)*wl[2];
			wr[3] = wl[3];
			break;
		case 1:
			// Outflow
			wr = wl;
			break;
		case 2:
			wr[0] = 2.4739;
			wr[1] = 2.5876;
			wr[2] = -0.5438;
			wr[3] = 2.2685;
			break;
		case 3:
			wr[0] = 1.4;
			wr[1] = 2.9;
			wr[2] = 0.0;
			wr[3] = 1.0;
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
	// Definisco il modello (con gamma=1.4)
	myModel model(1.4);
	// Definisco la mesh
	myMesh mesh;
	// Leggo la mesh
	Mesh::IO::MeshReader(mesh, meshfile);
	// Definisco il solutore per il mio modello
	mySolver solver(model, mesh);
	// Inizializzo alcuni parametri
	solver.setCFLmax(0.3);
	solver.setIC(init);
	solver.setBC(bc);
	// Inizializzo il solutore
	solver.init();
	solver.setDirectory("./data");
	// Passi temporali
	for (int i = 0; i <= 3000; ++i) {
		std::cout << "== Timestep " << i << " == currtime: " << std::setw(8) << solver.getCurrTime();
		std::cout << ", dt = " << std::setw(8) << solver.getCurrDt() << std::endl;
		solver.timestep();
		if (i%50 == 0) solver.framegrab(i/50, gnuplot, interpolated);
	}
	return 0;
}

