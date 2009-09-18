#ifndef _FINITEVOLUME_HPP
#define _FINITEVOLUME_HPP

// Libreria per la Mesh
#include <solvers/finitevolume/mesh_finitevolume_traits.hpp>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>

namespace ConservationLaw2D {
	/*! \namespace Solver
	\brief Namespace dei solutori */
	namespace Solver {
		using namespace std;
		/*! \class FiniteVolume
		\brief Solutore a Volumi Finiti per leggi di conservazione 2d */
		template <typename MODEL, typename NUMFLUX>
		class FiniteVolume {
		
			private:
				typedef typename MODEL::real_t							real_t;
				typedef typename MODEL::SolType							SolType;
				typedef typename Mesh::DefaultTraits<real_t,SolType>	Traits;
			public:
				// Mesh
				/*! \brief Mesh specifica per volumi finiti */
				typedef typename Traits::PolygonalMesh		FVMesh;
			private:
				// Puntatori
				typedef typename FVMesh::polygon_ptr		polygon_ptr;
				// Iteratori
				typedef typename FVMesh::polygon_it			p_it;
				typedef typename FVMesh::vertex_it			v_it;
				typedef typename FVMesh::hedge_it			e_it;
				// Circolatori
				typedef typename FVMesh::Vertex::PolygonCirculator	vp_cit;
				typedef typename FVMesh::Polygon::HEdgeCirculator	he_cit;
				typedef typename FVMesh::Polygon::VertexCirculator	pv_cit;
				// Condizioni iniziali e termine sorgente
				// [in ingresso x, y e colore]
				typedef SolType (*INITCOND)( size_t, real_t, real_t );
				// [in ingresso x, y, t e colore]
				typedef SolType (*BOUNDARYCOND)( SolType&, size_t, real_t, real_t, real_t, real_t, real_t );
				// [in ingresso x, y, t e colore]
				typedef SolType (*SOURCE)( const SolType&, size_t, real_t, real_t, real_t );
				
				// Valuta rhs del poligono dato
				inline SolType RHS( const polygon_ptr ) const;

			public:
				/*! \brief Costruttore del solutore 
				 \param[in] model Referenza al modello di problema
				 \param[in] mesh Referenza alla mesh
				*/
				FiniteVolume( MODEL& model, FVMesh& mesh )
					:model_(model),mesh_(mesh),NumFlux(model),cflmax_(0.0),hmax_(0.0),currtime_(0.0) {};
				
				// Impostazioni
				/*! \brief Imposta il massimo CFL */
				void setCFLmax ( real_t cflm ) { cflmax_ = cflm; }
				/*! \brief Imposta le condizioni iniziali */
				void setIC ( INITCOND ic ) { InitialCondition = ic; }
				/*! \brief Imposta le condizioni al bordo */
				void setBC ( BOUNDARYCOND bc ) { BoundaryCondition = bc; }
				/*! \brief Imposta il termine sorgente */
				void setSource ( SOURCE s ) { Source = s; }
				/*! \brief Imposta la directory nella quale sara' salvata la soluzione */
				void setDirectory ( const string& dir ) { datadir_ = dir; }
				// Accesso
				/*! \brief Restituisce il tempo corrente */
				real_t getCurrTime(void) { return currtime_; }
				/*! \brief Restituisce il passo temporale corrente */
				real_t getCurrDt(void) { return dt_; }
				
				// Inizializza il solutore
				/*! \brief Inizializza il solutore */
				void init();
				// Passo temporale
				/*! \brief Esegue un passo temporale */
				void timestep();
			private:
				void updateTimestep();
			public:
				/*! \brief Salva un frame della soluzione
				\param[in] id Id del frame
				\param[in] gnuplot Tipologia del plot (GNUPLOT, MATLAB o ALTRO)
				\param[in] interpolated Dati interpolati ai vertici */
				void framegrab(size_t const, bool, bool) const;

			private:
				// Modello
				MODEL& model_;
				// Mesh
				FVMesh& mesh_;
				// Flusso numerico
				NUMFLUX NumFlux;
				// Condizioni iniziali e bordo
				INITCOND		InitialCondition;
				BOUNDARYCOND	BoundaryCondition;
				SOURCE			Source;
				// Altro
				real_t cflmax_, hmax_, dt_, currtime_;
				string datadir_;
		};
		
		
		///////////////////////
		//  IMPLEMENTAZIONE  //
		///////////////////////
		
		template <typename MODEL,typename NUMFLUX>
		void FiniteVolume<MODEL,NUMFLUX>::init() {
			// Aggiorno Hmax e inizializzo la soluzione
			std::cout << "============================= " << std::endl;
			std::cout << "Init Finite Volume Solver ... " << std::endl;
			std::cout << "============================= " << std::endl;
			// Inizializzo la geometria per la mesh
			mesh_.init_geom();
			hmax_ = 0.0;
			for (p_it p = mesh_.p_begin(); p != mesh_.p_end(); ++p) {
				hmax_ = max( hmax_, (*p)->diam() );
				(*p)->sol = model_.PrimitiveToConservative(InitialCondition((*p)->getColor(),(*p)->cx(),(*p)->cy()));
			}
		}
		
		template <typename MODEL,typename NUMFLUX>
		inline typename MODEL::SolType FiniteVolume<MODEL,NUMFLUX>::RHS( const polygon_ptr p ) const {
			// Valuto il flusso attraverso i bordi
			SolType Flux = SolType::Zero();
			// Stato a sinistra
			SolType qlstate = p->sol0;
			// Stato a destra
			SolType qrstate;
			he_cit e = p->beginE();
			do {
				if ( !e->isBoundary() ) {
					// Lato interno
					qrstate = e->polygonR().sol0;
				} else {
					// Lato di bordo
					SolType wl = model_.ConservativeToPrimitive(qlstate);
					SolType wr = BoundaryCondition(wl,e->getColor(),e->xm(),e->ym(),e->nx(),e->ny(),currtime_);
					qrstate = model_.PrimitiveToConservative(wr);
				}
				Flux += e->length() * NumFlux(qlstate, qrstate, e->nx(), e->ny());
				++e;
			} while ( e != p->beginE() );
			// Valuto il termine sorgente, se presente
			SolType SourceTerm = SolType::Zero();
			//if ( !Source ) SourceTerm = Source( qlstate, p->getColor(), p->cx(), p->cy(), currtime_ );
			// Sommo e divido per l'area
			return (SourceTerm - Flux) / p->area();
		}

		template <typename MODEL, typename NUMFLUX>
		inline void FiniteVolume<MODEL,NUMFLUX>::timestep( void ) {
			// Salvo la soluzione al passo precedente e calcolo maxLambda
			updateTimestep();
			for (p_it p = mesh_.p_begin(); p != mesh_.p_end(); ++p) {
				(*p)->sol0 = (*p)->sol;
				if ( !model_.ConsistentState((*p)->sol0) ) {
					std::cerr << "Bad state solution! Maybe too high CFL number ..." << std::endl;
					exit(1);
				}
			}
			// Itero sui poligoni
			for (p_it p = mesh_.p_begin(); p != mesh_.p_end(); ++p) {
				// Risolvo l'ODE
				(*p)->sol = (*p)->sol0 + dt_ * RHS(*p);
			}
			// Aggiorno currtime_
			currtime_ += dt_;
		}
		
		template <typename MODEL, typename NUMFLUX>
		inline void FiniteVolume<MODEL,NUMFLUX>::updateTimestep(void) {
			// Aggiorno il passo temporale
			dt_ = 1e10;
			// Calcolo dt da CFL desiderato (cflmax)
			for (p_it p = mesh_.p_begin(); p != mesh_.p_end(); ++p) {
				dt_ = min( dt_, cflmax_* (*p)->diam()/model_.MaxLambda((*p)->sol) );
			}
		}
		
		template <typename MODEL, typename NUMFLUX>
		inline void FiniteVolume<MODEL,NUMFLUX>::framegrab( size_t const id, bool gnuplot, bool interpolated ) const {
			// Nome del file
			stringstream buffer;
			buffer.fill('0');
			buffer << datadir_ << "/solution" << std::setw(4) << id << ".dat";
			// Apro il file
			std::ofstream filehandle(buffer.str());
			// Scrivo i valori
			if (!gnuplot) {
				// Interpolo in ogni caso
				// Solo mesh triangolari
				assert(mesh_.isTriangular());
				for ( v_it v = mesh_.v_begin(); v != mesh_.v_end(); ++v ) {
					// Interpolo dai poligoni adiacenti
					SolType tmpsol = SolType::Zero();
					size_t count(0);
					vp_cit p = (*v)->beginP();
					do {
						tmpsol += model_.ConservativeToPrimitive(p->sol);
						count++;
						p++;
					} while( p != (*v)->beginP() );
					tmpsol /= count;
					// Output
					filehandle << tmpsol[0] << std::endl;
				}
			}
			if (gnuplot) {
				// Solo mesh triangolari
				assert(mesh_.isTriangular());
				// Itero su tutti i triangoli
				for ( p_it p = mesh_.p_begin(); p != mesh_.p_end(); ++p ) {
					if (!interpolated) {
						// Soluzione non interpolata
						SolType sol = model_.ConservativeToPrimitive((*p)->sol);
						stringstream strsol;
						for (int i=0; i<sol.rows(); ++i)
							strsol << sol[i] << " ";
						pv_cit vc = (*p)->beginV();
						// Vertice 1
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl;
						vc++;
						// Vertice 2
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl << endl;
						vc++;
						// Vertice 3
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl;
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl << endl << endl;
					} else {
						// Soluzione interpolata
						SolType tmpsol = SolType::Zero();
						size_t count;
						pv_cit vc;
						vp_cit pc;
						stringstream strsol;
						// Vertice 1
						vc = (*p)->beginV();
						count = 0.0;
						pc = vc->beginP();
						do {
							tmpsol += model_.ConservativeToPrimitive(pc->sol);
							count++;
							pc++;
						} while( pc != vc->beginP() );
						tmpsol /= count;
						for (int i=0; i<tmpsol.rows(); ++i)
							strsol << tmpsol[i] << " ";
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl;
						// Vertice 2
						vc++;
						tmpsol = SolType::Zero();
						count = 0.0;
						pc = vc->beginP();
						do {
							tmpsol += model_.ConservativeToPrimitive(pc->sol);
							count++;
							pc++;
						} while( pc != vc->beginP() );
						tmpsol /= count;
						strsol.str("");
						for (int i=0; i<tmpsol.rows(); ++i)
							strsol << tmpsol[i] << " ";
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl << endl;
						// Vertice 3
						vc++;
						tmpsol = SolType::Zero();
						count = 0.0;
						pc = vc->beginP();
						do {
							tmpsol += model_.ConservativeToPrimitive(pc->sol);
							count++;
							pc++;
						} while( pc != vc->beginP() );
						tmpsol /= count;
						strsol.str("");
						for (int i=0; i<tmpsol.rows(); ++i)
							strsol << tmpsol[i] << " ";
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl;
						filehandle << vc->x() << " " << vc->y() << " " << strsol.str() << endl << endl << endl;
					}
				}
			}
		}
	}
}

#endif
