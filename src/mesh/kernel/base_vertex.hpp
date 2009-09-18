#ifndef _MESH_BASE_VERTEX_HPP
#define _MESH_BASE_VERTEX_HPP

#include "base_kernel.hpp"

#include <vector>
#include <iostream>

#include "../iterators/circulatorTS.hpp"

namespace ConservationLaw2D {
	namespace Mesh {
		using std::vector;
		
		template <typename KERNEL>
		/*! \class BaseVertex
			\brief Definizione della struttura Vertex utilizzata nella Mesh.
		*/
		class BaseVertex : public KERNEL {
		
			typedef typename KERNEL::real_t			real_t;
			
			typedef typename KERNEL::Vertex			Vertex;
			typedef typename KERNEL::HEdge			HEdge;
			typedef typename KERNEL::Polygon		Polygon;
			
			typedef typename KERNEL::Vertex_ptr		vertex_ptr;
			typedef typename KERNEL::HEdge_ptr		hedge_ptr;
			typedef typename KERNEL::Polygon_ptr	polygon_ptr;
			
			friend class BaseHEdge<KERNEL>;
			friend class BasePolygon<KERNEL>;
			friend class BasePolygonalMesh<KERNEL>;
			
			public:
				// GEOMETRIA
				// Accesso alle coordinate
				/*! \brief Restituisce la coordinata \f$ x \f$ del vertice */
				inline real_t const & x(void) const { return x_; }
				/*! \brief Restituisce la coordinata \f$ x \f$ del vertice */
				inline real_t       & x(void)       { return x_; }
				/*! \brief Restituisce la coordinata \f$ y \f$ del vertice */
				inline real_t const & y(void) const { return y_; }
				/*! \brief Restituisce la coordinata \f$ y \f$ del vertice */
				inline real_t       & y(void)       { return y_; }
				// Imposta la posizione del vertice
				/*! \brief Imposta le coordinate \f$ (x,y) \f$ del vertice */
				void setPosition(const real_t x, const real_t y) { x_ = x; y_ = y; }
				
				// CIRCOLATORI
				/*! \brief Circolatore sui vertici adiacenti al vertice dato */
				typedef typename Circulators::CirculatorTS<KERNEL,Vertex,Vertex>	VertexCirculator;
				/*! \brief Circolatore sui lati adiacenti al vertice dato */
				typedef typename Circulators::CirculatorTS<KERNEL,HEdge,Vertex>		HEdgeCirculator;
				/*! \brief Circolatore sui poligoni adiacenti al vertice dato */
				typedef typename Circulators::CirculatorTS<KERNEL,Polygon,Vertex>	PolygonCirculator;
				
				/*! \brief Restituisce il circolatore sui vertici adiacenti */
				VertexCirculator	beginV(void) { return VertexCirculator( &hedges_ ); }
				/*! \brief Restituisce il circolatore sui lati adiacenti */
				HEdgeCirculator		beginE(void) { return HEdgeCirculator( &hedges_ ); }
				/*! \brief Restituisce il circolatore sui poligoni adiacenti */
				PolygonCirculator	beginP(void) { return PolygonCirculator( &hedges_ ); }
				
				// OUTPUT
				/*! \brief Output a video nella forma (x,y) */
				friend std::ostream & operator<<( std::ostream& OS, const Vertex& v ) {
					OS << "(" << v.x() << "," << v.y() << ")";
					return OS;
				}
			private:
				// DATA
				// Coordinate
				real_t x_, y_;
				// Lista HEdge per vertici degeneri
				vector<hedge_ptr> hedges_;
		};
	}
}

#endif
