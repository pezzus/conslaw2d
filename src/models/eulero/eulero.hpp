#ifndef EULERO_HPP
#define EULERO_HPP

// Modello equazioni di Eulero

#include <cstddef>
// Libreria EIGEN per l'Algebra
#include <Eigen/Core>
#include <cmath>

// Dimensione dello spazio di stato
// [rho, rho u, rho v, rho e]
#define DIMENSION 4

namespace ConservationLaw2D {
	/*! \namespace Model
	\brief Namespace che contiene le varie tipologie di modelli */
	namespace Model {
		
		template <typename T>
		/*! \class Eulero
		\brief Modello per le equazioni di Eulero 2D */
		class Eulero {
		
			public:
				// Defininzioni vettori
				// Vettore per la soluzione
				/*! \brief Tipo di dato reale, per esempio \c float o \c double */
				typedef T	real_t;
				/*! \brief Tipo per i vettori contenenti la soluzione */
				typedef Eigen::Matrix<T, DIMENSION, 1>	SolType;
				/*! \brief Tipo per le matrici contenenti il flusso \f$ m \times 2 \f$ */
				typedef Eigen::Matrix<T, DIMENSION, 2>	FluxType;
				/*! \brief Costruttore del modello */
				/*! \param[in] gamma Costante \f$ \gamma \f$ del gas */
				Eulero( real_t gamma ):GAMMA(gamma) {}
				
				// Accesso variabili
				/*! \brief Restituisce la densità \f$ \rho \f$ */
				inline real_t R( const SolType& q ) const { return q[0]; }
				/*! \brief Restituisce la prima componente della velocità, \f$ u \f$ */
				inline real_t U( const SolType& q ) const { return q[1]/q[0]; }
				/*! \brief Restituisce la seconda compontente della velocità, \f$ v \f$ */
				inline real_t V( const SolType& q ) const { return q[2]/q[0]; }
				/*! \brief Restituisce l'energia specifica \f$ e \f$ (cinetica più quella interna) */
				inline real_t E( const SolType& q ) const { return q[3]/q[0]; }
				/*! \brief Restituisce la pressione \f$ p \f$ */
				inline real_t P( const SolType& q ) const { return (GAMMA-1)*R(q)*( E(q)-0.5*(U(q)*U(q)+V(q)*V(q)) ); }
				/*! \brief Restituisce la velocità locale del suono \f$ c \f$ */
				inline real_t C( const SolType& q ) const { return sqrt( GAMMA*P(q)/R(q) ); }
				/*! \brief Restituisce l'entalpia \f$ h \f$ */
				inline real_t H( const SolType& q ) const { return E(q) + P(q)/R(q); }
				/*! \brief Restituisce la costante \f$ \gamma \f$ */
				inline real_t gamma(void) const { return GAMMA; }
				// Variabili primitive <-> Variabili conservate
				// [rho, u, v, p]      <-> [rho, rho u, rho v, rho e]
				/*! \brief Restituisce la soluzione nelle variabili conservate \f$ (\rho,\rho u, \rho v, \rho e) \f$  */
				inline SolType PrimitiveToConservative( const SolType& w ) const {
					SolType q = SolType::Zero();
					q[0] = w[0];
					q[1] = w[0]*w[1];
					q[2] = w[0]*w[2];
					q[3] = 0.5*w[0]*(w[1]*w[1] + w[2]*w[2]) + w[3]/(GAMMA-1);
					
					return q;
				}
				/*! \brief Restituisce la soluzione nelle variabili primitive \f$ (\rho, u, v, p) \f$  */
				inline SolType ConservativeToPrimitive( const SolType& q ) const {
					SolType w = SolType::Zero();
					w[0] = R(q);
					w[1] = U(q);
					w[2] = V(q);
					w[3] = P(q);
					
					return w;
				}
				
				// Controlla se lo stato e' consistente
				/*! \brief Controlla che lo stato del sistema sia consistente, cioè \f$ \rho > 0 \f$ e \f$ p > 0 \f$  */
				inline bool ConsistentState( const SolType& q ) const {
					// Pressione e densita' devono essere positivi!
					return ( P(q) > 0 && R(q) > 0 );
				}
				
				// Flusso
				/*! \brief Restituisce il flusso esatto  */
				inline FluxType Flux( const SolType& q ) {
					FluxType Flux = FluxType::Zero();
					// F(q)
					Flux(0,0) = q[1];
					Flux(1,0) = q[1]*q[1]/q[0] + P(q);
					Flux(2,0) = q[1]*q[2]/q[0];
					Flux(3,0) = ( q[3] + P(q) )*U(q);
					// G(q)
					Flux(0,1) = q[2];
					Flux(1,1) = q[1]*q[2]/q[0];
					Flux(2,1) = q[2]*q[2]/q[0] + P(q);
					Flux(3,1) = ( q[3] + P(q) )*V(q);
					
					return Flux;
				}
				
				// Flusso in direzione normale
				/*! \brief Restituisce il flusso esatto in direzione normale */
				inline SolType NormalFlux( const SolType& q, const real_t nx, const real_t ny ) const {
					SolType Flux_x = SolType::Zero();
					SolType Flux_y = SolType::Zero();
					
					Flux_x[0] = q[1];
					Flux_x[1] = q[1]*q[1]/q[0] + P(q);
					Flux_x[2] = q[1]*q[2]/q[0];
					Flux_x[3] = ( q[3] + P(q) )*U(q);
					
					Flux_y[0] = q[2];
					Flux_y[1] = q[1]*q[2]/q[0];
					Flux_y[2] = q[2]*q[2]/q[0] + P(q);
					Flux_y[3] = ( q[3] + P(q) )*V(q);
					
					return Flux_x * nx + Flux_y * ny;
				}
				
				// Massimo autovalore valutato in q
				/*! \brief Restituisce il massimo autovalore  */
				inline real_t MaxLambda( const SolType& q ) {
					return max(abs(U(q)),abs(V(q)))+C(q);
				}
				
				// Autovalori
				/*! \brief Restituisce gli autovalori del flusso esatto */
				inline SolType EigenValues( const SolType& q, const real_t nx, const real_t ny ) const {
					SolType Eig = SolType::Zero();

					real_t un =  U(q)*nx + V(q)*ny;

					Eig[0] = un - C(q);
					Eig[1] = un;
					Eig[2] = un;
					Eig[3] = un + C(q);
					
					return Eig;
				}
			private:
				real_t GAMMA;
		};
		
	}
}

#endif