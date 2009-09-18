#ifndef _MESH_BASE_HEDGE_HPP
#define _MESH_BASE_HEDGE_HPP

#include <iostream>
#include <cassert>
#include <cmath>

namespace ConservationLaw2D {
	namespace Mesh {
	
		template <typename KERNEL>
		/*!
			\class BaseHEdge
			\brief Definizione della struttura Half-Edge utilizzata nella Mesh.
		*/
		class BaseHEdge : public KERNEL {
			
			typedef typename KERNEL::real_t       real_t;
		
			typedef typename KERNEL::Vertex       Vertex;
			typedef typename KERNEL::HEdge        HEdge;
			typedef typename KERNEL::Polygon      Polygon;
			
			typedef typename KERNEL::Vertex_ptr   vertex_ptr;
			typedef typename KERNEL::HEdge_ptr    hedge_ptr;
			typedef typename KERNEL::Polygon_ptr  polygon_ptr;
			
			friend class BaseVertex<KERNEL>;
			friend class BasePolygon<KERNEL>;
			friend class BasePolygonalMesh<KERNEL>;

			public:
				// GEOMETRIA
				/*! \brief Calcola e restituisce la lunghezza del lato */
				inline real_t length(void) const { return sqrt(pow(vertexE().x()-vertexS().x(),2)+pow(vertexE().y()-vertexS().y(),2)); }
				/*! \brief Calcola e restituisce la coordinata \f$ x \f$ del punto medio */
				inline real_t xm(void) const { return 0.5 * ( vertexS().x() + vertexE().x() ); }
				/*! \brief Calcola e restituisce la coordinata \f$ y \f$ del punto medio */
				inline real_t ym(void) const { return 0.5 * ( vertexS().y() + vertexE().y() ); }
				/*! \brief Calcola e restituisce la prima componente della normale */
				inline real_t nx(void) const { return ( vertexE().y() - vertexS().y() ) / length(); }
				/*! \brief Calcola e restituisce la seconda componente della normale */
				inline real_t ny(void) const { return ( vertexS().x() - vertexE().x() ) / length(); }
				/*! \brief Calcola e restituisce la prima componente della tangente */
				inline real_t tx(void) const { return ( vertexE().x() - vertexS().x() ) / length(); }
				/*! \brief Calcola e restituisce la seconda componente della tangente */
				inline real_t ty(void) const { return ( vertexE().y() - vertexS().y() ) / length(); }
				
				// TOPOLOGIA
				/*! \brief Restituisce il vertice iniziale (secondo orientamento) */
				inline Vertex const  & vertexS(void) const      { return *vertex_; }
				/*! \brief Restituisce il vertice iniziale (secondo orientamento) */
				inline Vertex        & vertexS(void)            { return *vertex_; }
				/*! \brief Restituisce il vertice finale (secondo orientamento) */
				inline Vertex const  & vertexE(void) const      { return *(nexthedge_->vertex_); }
				/*! \brief Restituisce il vertice finale (secondo orientamento) */
				inline Vertex        & vertexE(void)            { return *(nexthedge_->vertex_); }
				/*! \brief Restituisce il poligono alla sua sinistra (secondo orientamento) */
				inline Polygon const & polygonL(void) const     { return *polygon_; }
				/*! \brief Restituisce il poligono alla sua sinistra (secondo orientamento) */
				inline Polygon       & polygonL(void)           { return *polygon_; }
				/*! \brief Restituisce il poligono alla sua destra (secondo orientamento)
				\warning Il poligono potrebbe non esserci se il lato è di bordo! */
				inline Polygon const & polygonR(void) const     { assert(!isBoundary()); return *(twinhedge_->polygon_); }
				/*! \brief Restituisce il poligono alla sua destra (secondo orientamento)  */
				inline Polygon       & polygonR(void)           { assert(!isBoundary()); return *(twinhedge_->polygon_); }
				/*! \brief Chiediamo se il lato è di bordo */
				inline bool isBoundary(void) const              { return (twinhedge_ == NULL); }
				/*! \brief Restituisce l'halfedge successivo */
				inline HEdge const   & getNextHEdge(void) const { return *nexthedge_; }
				/*! \brief Restituisce l'halfedge successivo */
				inline HEdge         & getNextHEdge(void)       { return *nexthedge_; }
				/*! \brief Restituisce l'halfedge gemello */
				inline HEdge const   & getTwinHEdge(void) const { return *twinhedge_; }
				/*! \brief Restituisce l'halfedge gemello */
				inline HEdge         & getTwinHEdge(void)       { return *twinhedge_; }
				/*! \brief Restituisce il colore del lato */
				size_t	getColor(void) 		const	{ return color_; }
				/*! \brief Assegna un colore al lato */
				void	setColor(size_t c) 			{ color_ = c; }
				/*! \brief Output a video nella forma [(x1,y1):(x2,y2)] */
				friend std::ostream & operator<<( std::ostream& OS, const HEdge& e ) {
					OS << "[" << e.vertexS() << ":" << e.vertexE() << "]";
					return OS;
				}
				
			private:
				// DATA
				// Vertice di partenza
				vertex_ptr	vertex_;
				// Poligono alla sinistra
				polygon_ptr	polygon_;
				// Halfedge successivo
				hedge_ptr	nexthedge_;
				// Halfedge gemello (puo' essere nullo)
				hedge_ptr	twinhedge_;
				// Colore
				size_t		color_;
		};
	}
}

#endif
