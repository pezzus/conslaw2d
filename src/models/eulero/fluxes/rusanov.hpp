#ifndef EULERO_RUSANOV_HPP
#define EULERO_RUSANOV_HPP

// Flusso approssimato di Rusanov

namespace ConservationLaw2D {
	namespace NumericalFlux {
		
		/*! \class Rusanov
		\brief Funtore per flusso numerico di Rusanov per modello di Eulero */
		template <typename MODEL>
		class Rusanov {
		
			typedef typename MODEL::real_t		real_t;
			typedef typename MODEL::SolType		SolType;
			typedef typename MODEL::FluxType	FluxType;
			public:
				/*! \brief Costruttore, prende in ingresso referenza al modello usato */
				Rusanov(MODEL& m):model(m) {}
				
				/*! \brief Resistuisce il flusso numerico
				\param[in] ql Stato del sistema a sinistra del lato
				\param[in] qr Stato del sistema a destra del lato 
				\param[in] nx Prima componente della normale al lato
				\param[in] ny Seconda componente della normale al lato 
				\return Flusso numerico di Rusanov nella direzione normale */
				inline SolType operator()(const SolType& ql, const SolType& qr, const real_t nx, const real_t ny) const {
					// Stato con momento in direzione normale e tangente
					SolType qql = ql;
					qql[1] =   ql[1]*nx + ql[2]*ny; // normale
					qql[2] = - ql[1]*ny + ql[2]*nx; // tangente
					SolType qqr = qr;
					qqr[1] =   qr[1]*nx + qr[2]*ny; // normale
					qqr[2] = - qr[1]*ny + qr[2]*nx; // tangente
					// Calcolo le variabili primitive nel sistema ruotato
					SolType wl = model.ConservativeToPrimitive(qql);
					SolType wr = model.ConservativeToPrimitive(qqr);
					// Calcolo i flussi nel sistema ruotato
					SolType FluxL = model.Flux(qql).col(0);
					SolType FluxR = model.Flux(qqr).col(0);
					// Calcolo le celerita'
					real_t cl = model.C(qql);
					real_t cr = model.C(qqr);
					// Calcolo S+
					real_t Splus = max( abs(wl[1])+cl, abs(wr[1])+cr );
					// NumericalFlux = (FR + FL)/2 - Splus*(qr-ql)/2
					SolType Flux = 0.5*( FluxL + FluxR ) - 0.5*Splus*( qqr - qql );
					// Torno alle variabili cartesiane
					SolType FFlux = Flux;
					FFlux[1] = Flux[1]*nx - Flux[2]*ny;
					FFlux[2] = Flux[1]*ny + Flux[2]*nx;
					return FFlux;
				}
			private:
				MODEL& model;
		};
	}
}

#endif