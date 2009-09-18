#ifndef SHALLOW_WATER_HPP
#define SHALLOW_WATER_HPP

// Modello equazioni Onde in acqua bassa

#include <cstddef>
// Libreria EIGEN per l'Algebra
#include <Eigen/Core>
#include <cmath>

// Dimensione dello spazio di stato
// [h, h*u, h*v]
#define DIMENSION 3
#define GRAVITY 9.81

namespace ConservationLaw2D {
	namespace Model {
		
		template <typename T>
		/*! \class ShallowWater
		\brief Modello per le equazioni di onde in acque basse */
		class ShallowWater {
		
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
				ShallowWater() {}
				
				// Flusso in direzione normale
				/*! \brief Restituisce il flusso esatto  */
				inline FluxType Flux( const SolType& q ) const {
					FluxType Flux = FluxType::Zero();
					// F(q)
					Flux(0,0) = q[1];
					Flux(1,0) = q[1]*q[1]/q[0] + 0.5*GRAVITY*q[0]*q[0];
					Flux(2,0) = q[1]*q[2]/q[0];
					// G(q)
					Flux(0,1) = q[2];
					Flux(1,1) = q[1]*q[2]/q[0];
					Flux(2,1) = q[2]*q[2]/q[0] + 0.5*GRAVITY*q[0]*q[0];

					return Flux;
				}
				
				// Variabili primitive <-> Variabili conservate
				/*! \brief Restituisce la soluzione nelle variabili conservate \f$ (\rho,\rho u, \rho v) \f$  */
				inline SolType PrimitiveToConservative( const SolType& w ) const {
					SolType q = SolType::Zero();
					q[0] = w[0];
					q[1] = w[0]*w[1];
					q[2] = w[0]*w[2];
					return q;
				}
				
				/*! \brief Restituisce la soluzione nelle variabili primitive \f$ (\rho, u, v) \f$  */
				inline SolType ConservativeToPrimitive( const SolType& q ) const {
					SolType w = SolType::Zero();
					w[0] = q[0];
					w[1] = q[1]/q[0];
					w[2] = q[2]/q[0];
					return w;
				}
				
				// Controlla se lo stato e' consistente
				/*! \brief Controlla che lo stato del sistema sia consistente, cioè se \f$ \rho > 0 \f$  */
				inline bool ConsistentState( const SolType& q ) const {
					// Densita' deve essere positiva!
					return ( q[0] > 0 );
				}
				
				// Massimo autovalore
				/*! \brief Restituisce il massimo autovalore  */
				inline real_t MaxLambda( const SolType& q ) {
					return pow(q[1]/q[0],2) + pow(q[2]/q[0],2) + sqrt(GRAVITY*q[0]);
				}
				
				// Autovalori
				/*! \brief Restituisce gli autovalori del flusso esatto */
				inline SolType EigenValues( const SolType& q, const real_t nx, const real_t ny ) const {
					SolType Eig = SolType::Zero();
					real_t U = q[1]/q[0];
					real_t V = q[2]/q[0];
					real_t C = sqrt(GRAVITY*q[0]);
					
					Eig[0] = U*nx + V*ny;
					Eig[1] = U*nx + V*ny - C;
					Eig[2] = U*nx + V*ny + C;

					return Eig;
				}
		};

	}
}

#endif
