#ifndef _CIRCULATOR_TS_IMPL_HPP
#define _CIRCULATOR_TS_IMPL_HPP

#include <cstddef>
#include <vector>

namespace ConservationLaw2D {
	namespace Mesh {
		namespace Circulators {

			////////////////////////
			// VERTEX CIRCULATORS //
			////////////////////////
			
			// VERTICE
			template<typename KERNEL>
			inline typename KERNEL::Vertex&
			getTS(const typename KERNEL::HEdge_ptr & current_hedge, 
			const typename KERNEL::Vertex& dummyT, const typename KERNEL::Vertex& dummyS)
			{ return current_hedge->getNextHEdge().vertexS(); }
			
			template <typename KERNEL>
			inline void nextTS(typename KERNEL::HEdge_ptr & current_hedge, const std::vector<typename KERNEL::HEdge_ptr>* hedges,
				std::size_t & current_vertex_hedge, bool & boundary, const typename KERNEL::Vertex& dummyT, 
				const typename KERNEL::Vertex& dummyS) {
				// Cerco il successivo
				typedef typename KERNEL::HEdge_ptr hedge_ptr;
				if (boundary) {
					// Caso lato degenere
					current_vertex_hedge = (current_vertex_hedge+1) % hedges->size();
					current_hedge = (*hedges)[current_vertex_hedge];
					boundary = false;
				} else {
					// Caso standard
					hedge_ptr tmphe(current_hedge);
					do {
						current_hedge = &(current_hedge->getNextHEdge());
					} while ( &(current_hedge->getNextHEdge().getNextHEdge()) != tmphe);
					if ( current_hedge->getNextHEdge().isBoundary() ) {
						// Halfedge di bordo
						boundary = true;
					} else {
						// Tutto regolare
						current_hedge = &(current_hedge->getNextHEdge().getTwinHEdge());
					}
				}
			}
			
			// POLIGONO
			template<typename KERNEL>
			inline typename KERNEL::Polygon&
			getTS(const typename KERNEL::HEdge_ptr & current_hedge, 
			const typename KERNEL::Polygon& dummyT, const typename KERNEL::Vertex& dummyS)
			{ return current_hedge->polygonL(); }

			template <typename KERNEL>
			inline void nextTS(typename KERNEL::HEdge_ptr & current_hedge, const std::vector<typename KERNEL::HEdge_ptr>* hedges,
				std::size_t & current_vertex_hedge, bool & boundary, 
				const typename KERNEL::Polygon& dummyT, const typename KERNEL::Vertex& dummyS) {
				// Cerco il successivo
				typedef typename KERNEL::HEdge_ptr hedge_ptr;
				// Procedura standard
				hedge_ptr tmphe(current_hedge);
				do {
					current_hedge = &(current_hedge->getNextHEdge());
				} while ( &(current_hedge->getNextHEdge()) != tmphe);
				if ( current_hedge->isBoundary() ) {
					// Halfedge di bordo
					current_vertex_hedge = (current_vertex_hedge+1) % hedges->size();
					current_hedge = (*hedges)[current_vertex_hedge];
				} else {
					// Tutto regolare
					current_hedge = &(current_hedge->getTwinHEdge());
				}
			}
			
			// HALF-EDGE
			template<typename KERNEL>
			inline typename KERNEL::HEdge&
			getTS(const typename KERNEL::HEdge_ptr & current_hedge, 
			const typename KERNEL::HEdge& dummyT, const typename KERNEL::Vertex& dummyS)
			{ return *current_hedge; }
			
			template <typename KERNEL>
			inline void nextTS(typename KERNEL::HEdge_ptr & current_hedge, const std::vector<typename KERNEL::HEdge_ptr>* hedges,
				std::size_t & current_vertex_hedge, bool & boundary, 
				const typename KERNEL::HEdge& dummyT, const typename KERNEL::Vertex& dummyS) {
				// Cerco il successivo
				typedef typename KERNEL::HEdge_ptr hedge_ptr;
				if (boundary) {
					// Caso lato degenere
					current_vertex_hedge = (current_vertex_hedge+1) % hedges->size();
					current_hedge = (*hedges)[current_vertex_hedge];
					boundary = false;
				} else {
					// Caso standard
					hedge_ptr tmphe(current_hedge);
					do {
						current_hedge = &(current_hedge->getNextHEdge());
					} while ( &(current_hedge->getNextHEdge()) != tmphe);
					if ( current_hedge->isBoundary() ) {
						// Halfedge di bordo
						boundary = true;
					} else {
						// Tutto regolare
						current_hedge = &(current_hedge->getTwinHEdge());
					}
				}
			}

			/////////////////////////
			// POLYGON CIRCULATORS //
			/////////////////////////
			
			// VERTICE
			template<typename KERNEL>
			inline typename KERNEL::Vertex&
			getTS(const typename KERNEL::HEdge_ptr & current_hedge, 
			const typename KERNEL::Vertex& dummyT, const typename KERNEL::Polygon& dummyS)
			{ return current_hedge->vertexS(); }
			
			template <typename KERNEL>
			inline void nextTS(typename KERNEL::HEdge_ptr & current_hedge, const std::vector<typename KERNEL::HEdge_ptr>* hedges,
				std::size_t & current_vertex_hedge, bool & boundary, const typename KERNEL::Vertex& dummyT, 
				const typename KERNEL::Polygon& dummyS) {
				// Semplicemente la successiva
				current_hedge = &(current_hedge->getNextHEdge());
			}
			
			// POLIGONO
			template<typename KERNEL>
			inline typename KERNEL::Polygon&
			getTS(const typename KERNEL::HEdge_ptr & current_hedge, 
			const typename KERNEL::Polygon& dummyT, const typename KERNEL::Polygon& dummyS)
			{ return current_hedge->polygonR(); }

			template <typename KERNEL>
			inline void nextTS(typename KERNEL::HEdge_ptr & current_hedge, const std::vector<typename KERNEL::HEdge_ptr>* hedges,
				std::size_t & current_vertex_hedge, bool & boundary, 
				const typename KERNEL::Polygon& dummyT, const typename KERNEL::Polygon& dummyS) {
				// Cerco il successivo
				typedef typename KERNEL::HEdge_ptr hedge_ptr;
				// Procedura standard: vado avanti finche' non trovo un lato non di bordo
				// Caso poligoni isolati: restituisce errore con assert alla prima lettura
				hedge_ptr tmphe(current_hedge);
				do {
					current_hedge = &(current_hedge->getNextHEdge());
				} while ( (current_hedge != tmphe) && current_hedge->isBoundary() );
			}
			
			// HALF-EDGE
			template<typename KERNEL>
			inline typename KERNEL::HEdge&
			getTS(const typename KERNEL::HEdge_ptr & current_hedge, 
			const typename KERNEL::HEdge& dummyT, const typename KERNEL::Polygon& dummyS)
			{ return *current_hedge; }
			
			template <typename KERNEL>
			inline void nextTS(typename KERNEL::HEdge_ptr & current_hedge, const std::vector<typename KERNEL::HEdge_ptr>* hedges,
				std::size_t & current_vertex_hedge, bool & boundary, 
				const typename KERNEL::HEdge& dummyT, const typename KERNEL::Polygon& dummyS) {
				// Semplicemente la successiva
				current_hedge = &(current_hedge->getNextHEdge());
			}
		}
	}
}

#endif