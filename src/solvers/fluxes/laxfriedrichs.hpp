#ifndef LAXFRIEDRICHS_HPP
#define LAXFRIEDRICHS_HPP

// Flusso Lax-Friedrichs
// Generico

namespace ConservationLaw2D {
	namespace NumericalFlux {
		
		// Lax-Friedrichs
		/*! \class LaxFriedrichs
		\brief Funtore per flusso numerico di Lax-Friedrichs generico */
		template <typename MODEL>
		class LaxFriedrichs {
		
			typedef typename MODEL::real_t		real_t;
			typedef typename MODEL::SolType		SolType;
			typedef typename MODEL::FluxType	FluxType;
			public:
				/*! \brief Costruttore, prende in ingresso referenza al modello usato */
				LaxFriedrichs(MODEL& m):model(m) {}
				
				/*! \brief Resistuisce il flusso numerico
				\param[in] ql Stato del sistema a sinistra del lato
				\param[in] qr Stato del sistema a destra del lato 
				\param[in] nx Prima componente della normale al lato
				\param[in] ny Seconda componente della normale al lato 
				\return Flusso numerico di Lax-Friedrichs nella direzione normale */
				inline SolType operator()(const SolType& ql, const SolType& qr, const real_t nx, const real_t ny) const { 
					// Calcolo i flussi
					FluxType FluxL = model.Flux(ql);
					FluxType FluxR = model.Flux(qr);
					// Viscosita' artificiale
					real_t sigma = max( model.MaxLambda(ql), model.MaxLambda(qr) );
					// Flusso in direzione normale
					return 0.5 * ( (FluxL.col(0)+FluxR.col(0))*nx + (FluxL.col(1)+FluxR.col(1))*ny - sigma*(qr-ql) );
				}
			private:
				MODEL& model;
		};
	}
}

#endif