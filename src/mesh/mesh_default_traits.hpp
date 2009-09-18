#ifndef _MESH_DEFAULT_TRAITS_HPP
#define _MESH_DEFAULT_TRAITS_HPP

#include <mesh/kernel/base_kernel.hpp>
#include <mesh/kernel/base_vertex.hpp>
#include <mesh/kernel/base_hedge.hpp>
#include <mesh/kernel/base_polygon.hpp>
#include <mesh/kernel/base_polygonalmesh.hpp>

namespace ConservationLaw2D {
	namespace Mesh {
		template <typename T>
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
			class HEdge : public BaseHEdge<Kernel> {};
			/*! \struct Polygon
				\brief Definizione di Polygon
			*/
			class Polygon : public BasePolygon<Kernel> {};
			/*! \struct PolygonalMesh
				\brief Mesh di tipo poligonale
			*/
			class PolygonalMesh : public BasePolygonalMesh<Kernel> {};
		};
	}
}

#endif
