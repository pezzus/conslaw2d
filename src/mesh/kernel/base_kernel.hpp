#ifndef _MESH_BASE_KERNEL_HPP
#define _MESH_BASE_KERNEL_HPP

namespace ConservationLaw2D {
	namespace Mesh {
		// Forward Declaration
		template <typename KERNEL> class BaseVertex;
		template <typename KERNEL> class BaseHEdge;
		template <typename KERNEL> class BasePolygon;
		template <typename KERNEL> class BasePolygonalMesh;

		template <typename RTYPE, typename VERTEX, typename HEDGE, typename POLYGON, typename MESH>
		/*! \struct BaseKernel
		\brief Contenitore di tipi comuni a tutti gli elementi della mesh
		*/
		struct BaseKernel {
			/*! \brief Tipo di dato reale, ad esempio \c real o \c double */
			typedef RTYPE     real_t;
			/*! \brief Vertex generico */
			typedef VERTEX    Vertex;
			/*! \brief Halfedge generico */
			typedef HEDGE     HEdge;
			/*! \brief Polygon generico */
			typedef POLYGON   Polygon;
			/*! \brief Mesh generica */
			typedef MESH      Mesh;
			/*! \brief Puntatore a tipo Vertex generico */
			typedef Vertex*   Vertex_ptr;
			/*! \brief Puntatore a tipo Halfedge generico */
			typedef HEdge*    HEdge_ptr;
			/*! \brief Puntatore a tipo Polygon generico */
			typedef Polygon*  Polygon_ptr;
		};
	}
}

#endif
