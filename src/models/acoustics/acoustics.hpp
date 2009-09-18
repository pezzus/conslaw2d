#ifndef ACOUSTICS_HPP
#define ACOUSTICS_HPP

// Modello equazioni Acustica

#include <cstddef>
// Libreria EIGEN per l'Algebra
#include <Eigen/Core>
#include <cmath>

// Dimensione dello spazio di stato
// [p, u, v]
#define DIMENSION 3


namespace ConservationLaw2D {
	namespace Model {
		
		/*! \class LinearAcoustics
		\brief Modello per le equazioni dell'acustica lineari */
		template <typename T>
		class LinearAcoustics {
		
			public:
				// Defininzioni vettori
				// Vettore per la soluzione
				/*! \brief Tipo di dato reale, per esempio \c float o \c double */
				typedef T	real_t;
				/*! \brief Tipo per i vettori contenenti la soluzione */
				typedef Eigen::Matrix<T, DIMENSION, 1>	SolType;
				/*! \brief Tipo per le matrici contenenti il flusso \f$ m \times 2 \f$ */
				typedef Eigen::Matrix<T, DIMENSION, 2>	FluxType;
				/*! \brief Costruttore del modello
				\param[in] rho0 Densità del gas nello stato non perturbato
				\param[in] K0 Modulo di elasticità lineare */
				LinearAcoustics( real_t rho0, real_t K0 )
					:rho0_(rho0), K0_(K0), c0_(sqrt(K0_/rho0_)), Z0_(c0_*rho0_) {}
				
				// Accesso ai dati
				/*! \brief Restituisce l'impedenza del mezzo \f$ Z_0 \f$ */
				real_t Z(void) { return Z0_; }
				/*! \brief Restituisce la velocità del suono nel mezzo \f$ c_0 \f$ */
				real_t C(void) { return c0_; }
				
				// Flusso in direzione normale
				/*! \brief Restituisce il flusso esatto  */
				inline FluxType Flux( const SolType& q ) const {
					FluxType Flux = FluxType::Zero();
					// F(q)
					Flux(0,0) = K0_ * q[1];
					Flux(1,0) = q[0]/rho0_;
					Flux(2,0) = 0.0;
					// G(q)
					Flux(0,1) = K0_ * q[2];
					Flux(1,1) = 0.0;
					Flux(2,1) = q[0]/rho0_;
					
					return Flux;
				}
				
				// Variabili primitive <-> Variabili conservate
				/*! \brief Restituisce la soluzione nelle variabili conservate \f$ (p,u,v) \f$  */
				inline SolType PrimitiveToConservative( const SolType& w ) const {
					return w;
				}
				
				/*! \brief Restituisce la soluzione nelle variabili primitive \f$ (p,u,v) \f$  */
				inline SolType ConservativeToPrimitive( const SolType& q ) const {
					return q;
				}
				
				// Controlla se lo stato e' consistente
				/*! \brief Controlla se lo stato è consistente  */
				inline bool ConsistentState( const SolType& q ) const {
					// Pressione positiva?
					return true;
				}
				
				// Massimo autovalore
				/*! \brief Restituisce il massimo autovalore  */
				inline real_t MaxLambda( const SolType& q ) {
					return c0_;
				}
				
				// Autovalori
				/*! \brief Restituisce gli autovalori del flusso esatto */
				inline SolType EigenValues( const SolType& q, const real_t nx, const real_t ny ) const {
					SolType Eig = SolType::Zero();

					Eig[0] = - c0_;
					Eig[1] = 0.0;
					Eig[2] = + c0_;
					
					return Eig;
				}
				
			private:
				real_t rho0_, K0_, c0_, Z0_;
		};

	}
}

#endif
