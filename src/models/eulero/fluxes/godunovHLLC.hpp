#ifndef EULERO_GODUNOV_HLLC_HPP
#define EULERO_GODUNOV_HLLC_HPP

// Flusso di Godunov
// Problema di Riemann risolto con approssimazione HLLC

namespace ConservationLaw2D {
	namespace NumericalFlux {
		
		/*! \class GodunovHLLC
		\brief Funtore per flusso numerico HLLC per il modello di Eulero */
		template <typename MODEL>
		class GodunovHLLC {
		
			typedef typename MODEL::real_t		real_t;
			typedef typename MODEL::SolType		SolType;
			typedef typename MODEL::FluxType	FluxType;
			public:
				/*! \brief Costruttore, prende in ingresso referenza al modello usato */
				GodunovHLLC(MODEL& m):model(m) {}
				
				/*! \brief Resistuisce il flusso numerico
				\param[in] ql Stato del sistema a sinistra del lato
				\param[in] qr Stato del sistema a destra del lato 
				\param[in] nx Prima componente della normale al lato
				\param[in] ny Seconda componente della normale al lato 
				\return Flusso numerico nella direzione normale, con approssimazione HLLC */
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
					// Zona Ustar
					real_t Sstar = (wr[3]-wl[3]+wl[0]*wl[1]*(SL-wl[1])-wr[0]*wr[1]*(SR-wr[1]))/(wl[0]*(SL-wl[1])-wr[0]*(SR-wr[1]));
					SolType qlstar( 1.0, Sstar, wl[2], model.E(qql)+(Sstar-wl[1])*(Sstar+wl[3]/(wl[0]*(SL-wl[1]))) );
					qlstar *= wl[0]*((SL-wl[1])/(SL-Sstar));
					SolType qrstar( 1.0, Sstar, wr[2], model.E(qqr)+(Sstar-wr[1])*(Sstar+wr[3]/(wr[0]*(SR-wr[1]))) );
					qrstar *= wr[0]*((SR-wr[1])/(SR-Sstar));
					SolType Flux;
					if ( SL >= 0 ) {
						Flux = FluxL;
					} else {
						if ( Sstar >= 0 ) {
							Flux = FluxL + SL*(qlstar-qql);
						} else {
							if ( SR >= 0 ) {
								Flux = FluxR + SR*(qrstar-qqr);
							} else {
								Flux = FluxR;
							}
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