#include <models/eulero/eulero.hpp>
#include <models/eulero/fluxes/godunovROE.hpp>
#include <solvers/finitevolume.hpp>
#include <mesh/io/meshreader.hpp>

#include <iostream>
#include <ctime>
#include <cmath>

#define GAMMA 1.4

using namespace std;
using namespace ConservationLaw2D;

typedef double real_t;
typedef Model::Eulero<real_t>                    myModel;
typedef NumericalFlux::GodunovRoe<myModel>       myNumFlux;
typedef Solver::FiniteVolume<myModel,myNumFlux>  mySolver;
typedef mySolver::FVMesh                         myMesh;

// Tipo per la soluzione (vettore d-dim)
typedef myModel::SolType	SolType;

// Stato iniziale in variabili (p,u,v)
inline SolType init( size_t color, real_t x, real_t y ) {
	SolType sol = SolType::Zero();
	// rho, u, v
	real_t rho_out = 1.0, p_out = 1.0, u_out = 0.0;
	real_t rho_in  = 0.1, p_in  = 1.0, u_in  = 0.0;
	real_t p_inf = 10.0;
	real_t s_shock = u_out + sqrt(GAMMA*p_out/rho_out)*sqrt((GAMMA+1)/(2*GAMMA)*(p_inf/p_out)+(GAMMA-1)/(2*GAMMA));
	real_t rho_inf = rho_out*((GAMMA+1)*pow(u_out-s_shock,2))/((GAMMA-1)*pow(u_out-s_shock,2)+2*GAMMA*p_out/rho_out);
	real_t u_inf = (1.0-rho_out/rho_inf)*s_shock + u_out*rho_out/rho_inf;
	sol[0] = ( x > 0.2 ) ? rho_out : rho_inf;
	sol[1] = ( x > 0.2 ) ? u_out : u_inf;
	sol[3] = ( x > 0.2 ) ? p_out : p_inf;
	// Bolla
	sol[0] = ( (x-0.5)*(x-0.5)+y*y < 0.2*0.2 ) ? rho_in : sol[0];
	sol[1] = ( (x-0.5)*(x-0.5)+y*y < 0.2*0.2 ) ? u_in : sol[1];
	sol[3] = ( (x-0.5)*(x-0.5)+y*y < 0.2*0.2 ) ? p_in : sol[3];
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
	// Definisco il modello
	myModel model(GAMMA);
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
	// Inizializzo il solutore
	solver.init();
	solver.setDirectory("./data");
	// Passi temporali
	for (int i = 0; i < 5000; ++i) {
		std::cout << "== Timestep " << i << " == currtime: " << std::setw(8) << solver.getCurrTime();
		std::cout << ", dt = " << std::setw(8) << solver.getCurrDt() << std::endl;
		solver.timestep();
		if (i%50 == 0) solver.framegrab(i/50, gnuplot, interpolated);
	}
	return 0;
}

