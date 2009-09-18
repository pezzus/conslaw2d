#ifndef EULERO_GODUNOV_HLL_HPP
#define EULERO_GODUNOV_HLL_HPP

// Flusso di Godunov
// Problema di Riemann risolto con approssimazione HLL

namespace ConservationLaw2D {
	namespace NumericalFlux {
		
		// GodunovHLL
		/*! \class GodunovHLL
		\brief Funtore per flusso numerico HLL per il modello di Eulero */
		template <typename MODEL>
		class GodunovHLL {
		
			typedef typename MODEL::real_t		real_t;
			typedef typename MODEL::SolType		SolType;
			typedef typename MODEL::FluxType	FluxType;
			public:
				/*! \brief Costruttore, prende in ingresso referenza al modello usato */
				GodunovHLL(MODEL& m):model(m) {}
				
				/*! \brief Resistuisce il flusso numerico
				\param[in] ql Stato del sistema a sinistra del lato
				\param[in] qr Stato del sistema a destra del lato 
				\param[in] nx Prima componente della normale al lato
				\param[in] ny Seconda componente della normale al lato 
				\return Flusso numerico nella direzione normale, con approssimazione HLL */
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
					// Calcolo l'entalpia
					real_t Hl = (model.E(qql)+model.P(qql)/wl[0]);
					real_t Hr = (model.E(qqr)+model.P(qqr)/wr[0]);
					// Calcolo delle medie di ROE
					real_t uROE = (sqrt(wr[0])*wr[1]+sqrt(wl[0])*wl[1])/(sqrt(wr[0])+sqrt(wl[0]));
					real_t vROE = (sqrt(wr[0])*wr[2]+sqrt(wl[0])*wl[2])/(sqrt(wr[0])+sqrt(wl[0]));
					real_t HROE = (sqrt(wr[0])*Hr+sqrt(wl[0])*Hl)/(sqrt(wr[0])+sqrt(wl[0]));
					real_t cROE = sqrt( (model.gamma()-1)*( HROE-0.5*(uROE*uROE+vROE*vROE)) );
					// Stime alla Davis-Einfeldt-Roe di SL e SR
					real_t SL = uROE-cROE;
					real_t SR = uROE+cROE;
					// Stime alla Davis
					//real_t SL = min(wl[1]-model.C(qql), wr[1]-model.C(qqr));
					//real_t SR = min(wl[1]+model.C(qql), wr[1]+model.C(qqr));
					SolType Flux;
					if ( SL >= 0 ) {
						Flux = FluxL;
					} else {
						if ( SR >= 0 ) {
							Flux = (SR*FluxL - SL*FluxR + SL*SR*(qqr-qql))/(SR-SL);
						} else {
							Flux = FluxR;
						}
					}
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