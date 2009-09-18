#ifndef _MESH_DEFAULT_TRAITS_HPP
#define _MESH_DEFAULT_TRAITS_HPP

#include <mesh/kernel/base_kernel.hpp>
#include <mesh/kernel/base_vertex.hpp>
#include <mesh/kernel/base_hedge.hpp>
#include <mesh/kernel/base_polygon.hpp>
#include <mesh/kernel/base_polygonalmesh.hpp>

#include <algorithm>

namespace ConservationLaw2D {
	namespace Mesh {
		template <typename T, typename SOLTYPE>
		/*! \struct DefaultTraits
			\brief In questa struttura vengono definiti le classi vere e proprie per la Mesh
		*/
		struct DefaultTraits {
			class Kernel;
			class Vertex;
			class HEdge;
			class Polygon;
			class PolygonalMesh;

			/*! \struct Kernel
				\brief Struttura comune a tutti gli elementi della mesh
			*/
			class Kernel : public BaseKernel<T, Vertex, HEdge, Polygon, PolygonalMesh> {};
			/*! \struct Vertex
				\brief Definizione di Vertex
			*/
			class Vertex : public BaseVertex<Kernel> {};
			/*! \struct HEdge
				\brief Definizione di HEdge
			*/
			class HEdge : public BaseHEdge<Kernel> {
				public:
					inline T length() { return length_; }
					inline T xm() { return xm_; }
					inline T ym() { return ym_; }
					inline T nx() { return nx_; }
					inline T ny() { return ny_; }
					// Inizializza la geometria
					/*! \brief Inizializzo le quantità geometriche dei lati */
					void init_geom(void) {
						length_ = BaseHEdge<Kernel>::length();
						xm_ = BaseHEdge<Kernel>::xm();
						ym_ = BaseHEdge<Kernel>::ym();
						nx_ = BaseHEdge<Kernel>::nx();
						ny_ = BaseHEdge<Kernel>::ny();
					}
				private:
					T length_, xm_, ym_, nx_, ny_;
			};
			/*! \struct Polygon
				\brief Definizione di Polygon
			*/
			class Polygon : public BasePolygon<Kernel> {
				public:
					EIGEN_MAKE_ALIGNED_OPERATOR_NEW
					/*! \brief Restituisce l'area del poligono */
					inline T area() { return area_; }
					/*! \brief Restituisce il diametro dell'elemento */
					inline T diam() { return diam_; }
					/*! \brief Restituisce l'ascissa del baricentro dell'elemento */
					inline T cx() { return cx_; }
					/*! \brief Restituisce l'ordinata del baricentro dell'elemento */
					inline T cy() { return cy_; }
					// Inizializza la geometria
					/*! \brief Inizializzo le quantità geometriche del poligono */
					void init_geom(void) {
						area_ = BasePolygon<Kernel>::area();
						diam_ = BasePolygon<Kernel>::diam();
						cx_ = BasePolygon<Kernel>::cx();
						cy_ = BasePolygon<Kernel>::cy();
					}
					/*! \brief Soluzione corrente nel poligono */
					SOLTYPE sol;
					/*! \brief Soluzione al passo precedente nel poligono */
					SOLTYPE sol0;
				private:
					T area_, diam_, cx_, cy_;
			};
			/*! \struct PolygonalMesh
				\brief Mesh di tipo poligonale
			*/
			class PolygonalMesh : public BasePolygonalMesh<Kernel> {
			
				typedef BasePolygonalMesh<Kernel> parent;
				
				public:
					// Inizializza la geometria per chiamate rapide
					/*! \brief Inizializzo le quantità geometriche degli elementi della mesh*/
					void init_geom(void) {
						// HEdge
						typedef typename parent::hedge_it heit;
						for ( heit i = parent::he_begin(); i != parent::he_end(); ++i) {
							(*i)->init_geom();
						}
						// Poligoni
						typedef typename parent::polygon_it polyit;
						for ( polyit i = parent::p_begin(); i != parent::p_end(); ++i) {
							(*i)->init_geom();
						}
					}
			};
		};
	}
}

#endif
