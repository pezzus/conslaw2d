#ifndef EULERO_GODUNOV_ROE_HPP
#define EULERO_GODUNOV_ROE_HPP

// Flusso di Godunov
// Problema di Riemann risolto con approssimazione di Roe

namespace ConservationLaw2D {
	namespace NumericalFlux {
		
		// Godunov-Roe
		/*! \class GodunovRoe
		\brief Funtore per flusso numerico Godunov con approssimazione Roe per il modello di Eulero */
		template <typename MODEL>
		class GodunovRoe {
		
			typedef typename MODEL::real_t		real_t;
			typedef typename MODEL::SolType		SolType;
			typedef typename MODEL::FluxType	FluxType;
			public:
				/*! \brief Costruttore, prende in ingresso referenza al modello usato */
				GodunovRoe(MODEL& m):model(m) {}
				
				/*! \brief Resistuisce il flusso numerico
				\param[in] ql Stato del sistema a sinistra del lato
				\param[in] qr Stato del sistema a destra del lato 
				\param[in] nx Prima componente della normale al lato
				\param[in] ny Seconda componente della normale al lato 
				\return Flusso numerico nella direzione normale, con approssimazione ROE */
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
					// Calcolo le medie di ROE
					real_t srhol = sqrt(qql[0]), srhor = sqrt(qqr[0]);
					real_t rhoM  = srhol*srhor;
					real_t uM    = ( srhol*wl[1]+srhor*wr[1] )/( srhol+srhor );
					real_t vM    = ( srhol*wl[2]+srhor*wr[2] )/( srhol+srhor );
					real_t HM    = ( srhol*Hl+srhor*Hr )/( srhol+srhor );
					real_t C2M   = (model.gamma()-1)*( HM - 0.5*( uM*uM + vM*vM ) );
					real_t CM    = sqrt(C2M);
					// Calcolo gli autovalori
					SolType lambdaM( uM-CM, uM, uM, uM+CM );
					// Calcolo degli autovettori
					SolType K0( 1, uM-CM, vM, HM-uM*CM );
					SolType K1( 1, uM, vM, 0.5*(uM*uM+vM*vM) );
					SolType K2( 0, 0, 1, vM );
					SolType K3( 1, uM+CM, vM, HM+uM*CM );
					// Calcolo degli alpha
					SolType alphaM;
					alphaM[0] = 0.5/C2M*((wr[3]-wl[3])-rhoM*CM*(wr[1]-wl[1]));
					alphaM[1] = (wr[0]-wl[0]) - (wr[3]-wl[3])/C2M;
					alphaM[2] = rhoM*(wr[2]-wl[2]);
					alphaM[3] = 0.5/C2M*((wr[3]-wl[3])+rhoM*CM*(wr[1]-wl[1]));
					// Calcolo del flusso
					SolType Flux = SolType::Zero();
					if ( lambdaM[0] >= 0 ) {
						// Tutti positivi
						Flux = FluxL;
					} else {
						if ( lambdaM[1] >= 0 ) {
							// lambda1 < 0
							Flux = FluxL + alphaM[0]*lambdaM[0]*K0;
						} else {
							if ( lambdaM[3] >= 0 ) {
								// lambda1, lambda23 < 0, lambda4 > 0
								Flux = FluxR - alphaM[3]*lambdaM[3]*K3;
							} else {
								Flux = FluxR;
							}
						}
					}
					// Entropy fix
					// Calcolo dello stato nella regione star per entropy fix
					// Metodo ROE-Averaged
					real_t rhoLStar = wl[0] + alphaM[0];
					real_t uStar = (wl[0]*wl[1] + alphaM[0]*(uM-CM))/(wl[0]+alphaM[0]);
					real_t pStar = (model.gamma()-1)*(qql[3]+alphaM[0]*(HM-uM*CM)-0.5*rhoLStar*uStar*uStar);
					real_t CLStar = sqrt( model.gamma()*pStar/rhoLStar );
					real_t lambda1L = wl[1]-model.C(qql);
					real_t lambda1R = uStar-CLStar;
					if ( (lambda1L < 0) && (lambda1R > 0) ) {
						// Left transonic rarefaction wave
						Flux = FluxL + lambda1L*((lambda1R-lambdaM[0])/(lambda1R-lambda1L))*alphaM[0]*K0;
					}
					real_t rhoRStar = wr[0] - alphaM[3];
					uStar = (wr[0]*wr[1] + alphaM[3]*(uM+CM))/(wr[0]-alphaM[3]);
					pStar = (model.gamma()-1)*(qqr[3]-alphaM[3]*(HM+uM*CM)-0.5*rhoRStar*uStar*uStar);
					real_t CRStar = sqrt( model.gamma()*pStar/rhoRStar );
					real_t lambda4L = uStar+CRStar;
					real_t lambda4R = wr[1]+model.C(qqr);
					if ( (lambda4L < 0) && (lambda4R > 0) ) {
						// Right transonic rarefaction wave
						Flux = FluxR - lambda4R*((lambdaM[3]-lambda4L)/(lambda4R-lambda4L))*alphaM[3]*K3;
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