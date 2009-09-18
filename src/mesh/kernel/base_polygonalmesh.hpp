#ifndef _MESH_BASE_POLYGONALMESH_HPP
#define _MESH_BASE_POLYGONALMESH_HPP

#include <vector>
#include <iostream>
#include <algorithm>

namespace ConservationLaw2D {
	namespace Mesh {
		
		using namespace std;

		template <typename KERNEL>
		/*! \class BasePolygonalMesh
			\brief Definizione della classe primitiva della Mesh
		*/
		class BasePolygonalMesh : public KERNEL {
			public:
				/*! \brief Tipo di dato per numeri reali, per esempio \c double o \c float */
				typedef typename KERNEL::real_t			real_t;

				/*! \brief Tipo Vertex */
				typedef typename KERNEL::Vertex			Vertex;
				/*! \brief Tipo HEdge */
				typedef typename KERNEL::HEdge			HEdge;
				/*! \brief Tipo Polygon */
				typedef typename KERNEL::Polygon		Polygon;
			
				/*! \brief Puntatore a tipo Vertex */
				typedef typename KERNEL::Vertex_ptr		vertex_ptr;
				/*! \brief Puntatore a tipo HEdge */
				typedef typename KERNEL::HEdge_ptr		hedge_ptr;
				/*! \brief Puntatore a tipo Polygon */
				typedef typename KERNEL::Polygon_ptr	polygon_ptr;
			
			private:
				friend class BaseVertex<KERNEL>;
				friend class BaseHEdge<KERNEL>;
				friend class BasePolygon<KERNEL>;
			
				// Liste
				typedef vector<vertex_ptr> vertex_list;
				typedef vector<hedge_ptr> hedge_list;
				typedef vector<polygon_ptr> polygon_list;
				
				// Eliminazione puntatori
				template <typename T>
				static bool deleteAll( T* elm ) { delete elm; return true; }

			public:
				// Iteratori
				/*! \brief Iteratore sui vertici */
				typedef typename vertex_list::iterator		vertex_it;
				/*! \brief Iteratore sui lati */
				typedef typename hedge_list::iterator		hedge_it;
				/*! \brief Iteratore sui poligoni */
				typedef typename polygon_list::iterator		polygon_it;
			
			public:
				// ========================
				// CONSTRUCTOR & DESTRUCTOR
				/*! \brief Costruttore della mesh */
				BasePolygonalMesh():isTriangular_(true) {};
				/*! \brief Distruttore della mesh, si preoccupa di liberare la memoria */
				~BasePolygonalMesh() {
					// Elimino tutto
					std::remove_if(vertices_.begin(), vertices_.end(), deleteAll<Vertex>);
					std::remove_if(hedges_.begin(), hedges_.end(), deleteAll<HEdge>);
					std::remove_if(polygons_.begin(), polygons_.end(), deleteAll<Polygon>);
				};
				
				// MODIFICA MESH
				/*! \brief Metodo per l'aggiunta di un vertice alla mesh */
				inline vertex_ptr	addVertex	( const real_t x, const real_t y );
				/*! \brief Metodo per l'aggiunta di un poligono alla mesh
				\param[in] v Lista dei puntatori ai vertici */
				inline polygon_ptr	addPolygon	( const std::vector<vertex_ptr>& );
				
				// LETTURA MESH
				/*! \brief Restituisce il numero dei vertici */
				size_t nV(void) const { return vertices_.size(); }
				/*! \brief Restituisce il numero degli half-edge */
				size_t nE(void) const { return hedges_.size(); }
				/*! \brief Restituisce il numero dei poligoni */
				size_t nP(void) const { return polygons_.size(); }
				
				/*! \brief Restituisce l'i-esimo vertice */
				vertex_ptr v(size_t i) const { return vertices_[i]; }
				/*! \brief Restituisce l'i-esimo half-edge */
				hedge_ptr e(size_t i) const { return hedges_[i]; }
				/*! \brief Restituisce l'i-esimo poligono */
				polygon_ptr p(size_t i) const { return polygons_[i]; }

				/*! \brief Restituisce l'iteratore al primo vertice */
				vertex_it v_begin( void ) { return vertices_.begin(); }
				/*! \brief Restituisce l'iteratore al primo vertice */
				vertex_it v_begin( void ) const { return vertices_.begin(); }
				/*! \brief Restituisce l'iteratore all'ultimo vertice */
				vertex_it v_end( void ) { return vertices_.end(); }
				/*! \brief Restituisce l'iteratore all'ultimo vertice */
				vertex_it v_end( void ) const { return vertices_.end(); }
				
				/*! \brief Restituisce l'iteratore al primo half-edge */
				hedge_it he_begin( void ) { return hedges_.begin(); }
				/*! \brief Restituisce l'iteratore al primo half-edge */
				hedge_it he_begin( void ) const { return hedges_.begin(); }
				/*! \brief Restituisce l'iteratore all'ultimo half-edge */
				hedge_it he_end( void ) { return hedges_.end(); }
				/*! \brief Restituisce l'iteratore all'ultimo half-edge */
				hedge_it he_end( void ) const { return hedges_.end(); }
				
				/*! \brief Restituisce l'iteratore al primo poligono */
				polygon_it p_begin( void ) { return polygons_.begin(); }
				/*! \brief Restituisce l'iteratore al primo poligono */
				polygon_it p_begin( void ) const { return polygons_.begin(); }
				/*! \brief Restituisce l'iteratore all'ultimo poligono */
				polygon_it p_end( void ) { return polygons_.end(); }
				/*! \brief Restituisce l'iteratore all'ultimo poligono */
				polygon_it p_end( void ) const { return polygons_.end(); }
				
				/*! \brief Output di alcune statistiche sulla mesh */
				void stats( void );
				
				/*! \brief Restituisce vero se la mesh contiene solo triangolli */
				bool isTriangular() { return isTriangular_; }
				
			private:
				// DATA
				vertex_list		vertices_;
				hedge_list		hedges_;
				polygon_list	polygons_;
				bool			isTriangular_;
		};
		
		
		// ===========
		// DEFINITIONS
		// ===========
		
		template <typename KERNEL>
		inline typename KERNEL::Vertex_ptr BasePolygonalMesh<KERNEL>::addVertex(const real_t x, const real_t y) {
			// Creo il vertice nella memoria
			vertex_ptr vhandle( new Vertex() );
			// Imposto le coordinate
			vhandle->setPosition(x, y);
			// Salvo nella lista dei vertici
			vertices_.push_back( vhandle );
			// Restituisco l'handle
			return vhandle;
		}

		template <typename KERNEL>
		inline typename KERNEL::Polygon_ptr BasePolygonalMesh<KERNEL>::addPolygon( const std::vector<vertex_ptr>& v ) {
			// Creo il poligono
			polygon_ptr poly( new Polygon() );
			// Creo i lati del poligono
			size_t nsides = v.size();
			isTriangular_ &= (nsides == 3);
			vector<hedge_ptr> he(nsides);
			for (size_t i = 0; i < nsides; ++i) {
				// Nuovo halfhedge
				he[i] = (hedge_ptr)new HEdge();
				// Imposto il poligono a sinistra del lato
				he[i]->polygon_ = poly;
				// Imposto il vertice iniziale
				he[i]->vertex_ = v[i];
				// Aggiungo alla lista dei lati
				hedges_.push_back(he[i]);
			}
			// Aggiungo al poligono uno dei suoi lati
			poly->hedge_ = he[0];
			// Collegamenti tra halfedge
			typedef typename Vertex::HEdgeCirculator HECirc;
			for (size_t i = 0; i < nsides; ++i) {
				// Imposto il successivo
				he[i]->nexthedge_ = he[(i+1)%nsides];
				// Verifico se il vertice corrispondente e' nuovo
				if ( v[i]->hedges_.size() != 0 ) {
					// Gia' presente, cerco eventuali collegamenti
					HECirc vhec = v[i]->beginE();
					bool found = false;
					do {
						++vhec;
						found = ( vhec->vertex_ == v[(i+1)%nsides]);
					} while(!found && vhec != v[i]->beginE());
					if (found) {
						// Trovata corrispondenza, salvo
						he[i]->twinhedge_ = &(*vhec);
						poly->hedge_ = he[i];
						vhec->polygon_->hedge_ = &(*vhec);
					}
				}
			}
			// Collego il nuovo elemento con la mesh
			for (size_t i = 0; i < nsides; ++i) {
				v[i]->hedges_.push_back(he[i]);
				if ( he[i]->twinhedge_ ) {
					he[i]->twinhedge_->twinhedge_ = he[i];
				}
			}
			// Aggiorno la struttura dati dei vertici
			typedef typename vector<hedge_ptr>::iterator Iterator;
			for (size_t i = 0; i < nsides; ++i) {
				Iterator it = v[i]->hedges_.begin();
				while ( it != v[i]->hedges_.end() ) {
					if ( !(*it)->isBoundary() ) {
						it = v[i]->hedges_.erase( it );
					} else {
						++it;
					}
				}
				// Se e' vuoto ne metto uno a caso
				if ( v[i]->hedges_.empty() ) {
					v[i]->hedges_.push_back(he[i]);
				}
			}
			// Aggiungo il poligono alla mesh
			polygons_.push_back(poly);
			// Resituisco il puntatore
			return poly;
		}
		
		template <typename KERNEL>
		void BasePolygonalMesh<KERNEL>::stats( void ) {
			cout << "Mesh stats:" << endl;
			cout << " Number of vertices: " << nV() << endl;
			cout << " Number of hedges: " << nE() << endl;
			cout << " Number of polygons: " << nP() << endl;
			cout << " isTriangular: " << ((isTriangular_)?"yes":"no") << endl;
			cout << "Memory stats:" << endl;
			cout << " Vertices: " << (nV() * sizeof(Vertex))/1024. << " Kbytes" << endl;
			cout << " Edges: " << (nE() * sizeof(HEdge))/1024. << " Kbytes" << endl;
			cout << " Vertices: " << (nP() * sizeof(Polygon))/1024. << " Kbytes" << endl;
		}
	}
}

#endif
