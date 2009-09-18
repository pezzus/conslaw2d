#ifndef _MESH_BASE_POLYGON_HPP
#define _MESH_BASE_POLYGON_HPP

#include <cmath>

#include "../iterators/circulatorTS.hpp"

namespace ConservationLaw2D {
	namespace Mesh {
	
		template <typename KERNEL>
		/*! \class BasePolygon
			\brief Definizione della struttura Polygon utilizzata nella Mesh.
		*/
		class BasePolygon : public KERNEL {
			
			typedef typename KERNEL::real_t			real_t;
		
			typedef typename KERNEL::Vertex			Vertex;
			typedef typename KERNEL::HEdge			HEdge;
			typedef typename KERNEL::Polygon		Polygon;
			
			typedef typename KERNEL::Vertex_ptr		vertex_ptr;
			typedef typename KERNEL::HEdge_ptr		hedge_ptr;
			typedef typename KERNEL::Polygon_ptr	polygon_ptr;
			
			friend class BaseVertex<KERNEL>;
			friend class BaseHEdge<KERNEL>;
			friend class BasePolygonalMesh<KERNEL>;

			public:
				// GEOMETRIA
				/*! \brief Calcola e restituisce l'area del poligono */
				inline real_t area(void)	const;
				/*! \brief Calcola e restituisce l'ascissa del baricentro del poligono */
				inline real_t cx(void)		const;
				/*! \brief Calcola e restituisce l'ordinata del baricentro del poligono */
				inline real_t cy(void)		const;
				/*! \brief Calcola e restituisce il diametro del poligono */
				inline real_t diam(void)	const;
				
				// TOPOLOGIA
				/*! \brief Restituisce il colore del poligono */
				size_t	getColor(void) 		const	{ return color_; }
				/*! \brief Imposta il colore del poligono */
				void	setColor(size_t c) 			{ color_ = c; }

				// CIRCOLATORI
				/*! \brief Circolatore sui vertici adiacenti al poligono dato */
				typedef typename Circulators::CirculatorTS<KERNEL,Vertex,Polygon>	VertexCirculator;
				/*! \brief Circolatore sui lati adiacenti al poligono dato */
				typedef typename Circulators::CirculatorTS<KERNEL,HEdge,Polygon>	HEdgeCirculator;
				/*! \brief Circolatore sui poligoni adiacenti al poligono dato */
				typedef typename Circulators::CirculatorTS<KERNEL,Polygon,Polygon>	PolygonCirculator;
				
				/*! \brief Restituisce il circolatore sui vertici adiacenti */
				VertexCirculator	beginV(void) { return VertexCirculator( hedge_ ); }
				/*! \brief Restituisce il circolatore sui lati adiacenti */
				HEdgeCirculator		beginE(void) { return HEdgeCirculator( hedge_ ); }
				/*! \brief Restituisce il circolatore sui poligoni adiacenti */
				PolygonCirculator	beginP(void) { return PolygonCirculator( hedge_ ); }

				// OUTPUT
				/*! \brief Output a video nella forma {[lato_1] [lato_2] ... [lato_n]} */
				friend std::ostream & operator<<( std::ostream& OS, const Polygon& p ) {
					hedge_ptr tmphe = p.hedge_;
					OS << "{ ";
					do {
						OS << tmphe->vertexS() << " ";
						tmphe = &(tmphe->getNextHEdge());
					} while ( tmphe != p.hedge_ );
					OS << "}";
					return OS;
				}

			private:
				// DATA
				hedge_ptr	hedge_;
				size_t		color_;
		};
		
		// IMPLEMENTAZIONE //
		template <typename KERNEL>
		inline typename KERNEL::real_t BasePolygon<KERNEL>::area(void) const {
			// Calcolo dell'area di un poligono generico
			real_t area(0);
			hedge_ptr hei = hedge_;
			do {
				area += (hei->vertexS().x())*(hei->vertexE().y())-(hei->vertexE().x())*(hei->vertexS().y());
				hei = &(hei->getNextHEdge());
			} while ( hei != hedge_ );
			return area/2.0;
		}

		template <typename KERNEL>
		inline typename KERNEL::real_t BasePolygon<KERNEL>::cx(void) const {
			// Calcolo ascissa x del baricentro del poligono
			real_t cx(0.0);
			size_t nE(0);
			hedge_ptr hei = hedge_;
			do {
				cx += hei->vertexS().x();
				nE++;
				hei = &(hei->getNextHEdge());
			} while ( hei != hedge_ );
			return cx/nE;
		}

		template <typename KERNEL>
		inline typename KERNEL::real_t BasePolygon<KERNEL>::cy(void) const {
			// Calcolo ascissa y del baricentro del poligono
			real_t cy(0.0);
			size_t nE(0);
			hedge_ptr hei = hedge_;
			do {
				cy += hei->vertexS().y();
				nE++;
				hei = &(hei->getNextHEdge());
			} while ( hei != hedge_ );
			return cy/nE;
		}

		template <typename KERNEL>
		inline typename KERNEL::real_t BasePolygon<KERNEL>::diam(void) const {
			// Calcolo il diametro del poligono
			// Calcolo gli (n^2-n)/2 modi per unire i vertici e poi prendo il max
			real_t diam(0.0);
			hedge_ptr hei = hedge_;
			do {
				hedge_ptr hej = &(hei->getNextHEdge());
				do {
					diam = max(diam, pow(hei->vertexS().x()-hej->vertexS().x(),2)+pow(hei->vertexS().y()-hej->vertexS().y(),2));
					hej = &(hej->getNextHEdge());
				} while ( hej != hedge_ );
				hei = &(hei->getNextHEdge());
			} while ( hei != hedge_ );
			return sqrt(diam);
		}

	}
}

#endif
