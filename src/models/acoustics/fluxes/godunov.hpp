#ifndef ACOUSTICS_GODUNOV_HPP
#define ACOUSTICS_GODUNOV_HPP

// Flusso di Godunov
// Problema di Riemann risolto esattamente

namespace ConservationLaw2D {
	namespace NumericalFlux {
		
		// Godunov-Exact
		/*! \class Godunov
		\brief Funtore per flusso numerico di Godunov esatto per Acustica */
		template <typename MODEL>
		class Godunov {
		
			typedef typename MODEL::real_t		real_t;
			typedef typename MODEL::SolType		SolType;
			typedef typename MODEL::FluxType	FluxType;
			
			public:
				/*! \brief Costruttore, prende in ingresso referenza al modello usato */
				Godunov(MODEL& m):model(m) {}
				
				/*! \brief Resistuisce il flusso numerico
				\param[in] ql Stato del sistema a sinistra del lato
				\param[in] qr Stato del sistema a destra del lato 
				\param[in] nx Prima componente della normale al lato
				\param[in] ny Seconda componente della normale al lato 
				\return Flusso numerico di Godunov esatto nella direzione normale */
				inline SolType operator()(const SolType& ql, const SolType& qr, const real_t nx, const real_t ny) const {
					// Risolvo esattamente il problema di Riemann valutato in x/t = 0
					SolType delta = qr - ql;
					SolType r1( -model.Z(), nx, ny);
					real_t alpha = 0.5 * ( - 1.0/model.Z() * delta[0] + nx * delta[1] + ny * delta[2] );
					// Calcolo il flusso numerico
					FluxType Flux = model.Flux( ql+alpha*r1 );
					return Flux.col(0)*nx + Flux.col(1)*ny;
				}
			private:
				MODEL& model;
		};
	}
}

#endif