#ifndef _CIRCULATOR_TS_HPP
#define _CIRCULATOR_TS_HPP

#include <iterator>
#include <vector>

#include "circulatorTS_impl.hpp"

/*! \namespace ConservationLaw2D
\brief Namespace contenente tutte le classi e metodi per risolvere leggi di conservazione
*/
namespace ConservationLaw2D {
	/*! \namespace Mesh
	\brief Namespace contenente tutte le classi e metodi per gestire la mesh
	*/
	namespace Mesh {
		/*! \namespace Circulators
		\brief Namespace per i circolatori
		*/
		namespace Circulators {
		
			using std::vector;
			
			/*! \class CirculatorTS
				\brief Definizione del circolatore generico sull'oggetto T rispetto ad S
			*/
			template <typename KERNEL, typename T, typename S>
			class CirculatorTS : 
				public KERNEL,
				public std::iterator< std::forward_iterator_tag, T, size_t, T*, T& > 
			{
				private:
					typedef typename KERNEL::Vertex		Vertex;
					typedef typename KERNEL::HEdge		HEdge;
					typedef typename KERNEL::Polygon	Polygon;
			
					typedef typename KERNEL::Vertex_ptr		vertex_ptr;
					typedef typename KERNEL::HEdge_ptr		hedge_ptr;
					typedef typename KERNEL::Polygon_ptr	polygon_ptr;
					
					typedef typename std::iterator< std::bidirectional_iterator_tag, T, std::size_t, T*, T& > parent;
					typedef typename parent::value_type		value_type;
					typedef typename parent::reference		reference;
					typedef typename parent::pointer		pointer;
				public:
					/*! \brief Costruttore di default */
					CirculatorTS( void );
					/*! \brief Costruttore a partire da una lista di half-edge
					\warning Questo costruttore è utilizzato per gestire i vertici degeneri. */
					CirculatorTS( vector<hedge_ptr>* );
					/*! \brief Costruttore a partire da un half-edge singolo
					\warning Questo costruttore è utilizzato per gestire vertici regolari. */
					CirculatorTS( hedge_ptr );
					/*! \brief Costruttore copia */
					CirculatorTS( const CirculatorTS<KERNEL,T,S>& );
					/*! \brief Copia di un circolatore */
					CirculatorTS<KERNEL,T,S>& operator=( const CirculatorTS<KERNEL,T,S>& );
					
					/*! \brief Restituisce una referenza all'oggetto corrente */
					reference 	operator* ( void );
					/*! \brief Restituisce un puntatore all'oggetto corrente */
					pointer 	operator->( void );
					
					/*! \brief Passa all'oggetto successivo della lista */
					CirculatorTS<KERNEL,T,S>& operator++( void );
					/*! \brief Passa all'oggetto successivo della lista */
					CirculatorTS<KERNEL,T,S>  operator++( int );

					/*! \brief Controlla se i due circolatori puntano alla stessa cosa */
					bool operator==( const CirculatorTS<KERNEL,T,S>& rhs ) { return (current_hedge_ == rhs.current_hedge_); }
					/*! \brief Controlla se i due circolatori non puntano alla stessa cosa */
					bool operator!=( const CirculatorTS<KERNEL,T,S>& rhs ) { return (current_hedge_ != rhs.current_hedge_); }

				private:
					vector<hedge_ptr>*	hedges_;
					bool				boundary_;
					size_t				current_vertex_hedge_id_;
					hedge_ptr			current_hedge_;
					value_type			dummyT;
					S					dummyS;
			};
			
			// DEFINIZIONI
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S>::CirculatorTS( void ) 
				: hedges_(0), boundary_(false), current_vertex_hedge_id_(0), current_hedge_(0) {}
			
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S>::CirculatorTS( vector<typename KERNEL::HEdge_ptr>* ev )
				: hedges_(ev), boundary_(false)
			{
				current_vertex_hedge_id_ = (*hedges_).size()-1;
				current_hedge_ = (*hedges_)[current_vertex_hedge_id_];
			}
			
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S>::CirculatorTS( typename KERNEL::HEdge_ptr e )
				: hedges_(0), boundary_(false), current_vertex_hedge_id_(0), current_hedge_(e) { }
			
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S>::CirculatorTS( const CirculatorTS<KERNEL,T,S>& other ) 
				: hedges_(other.hedges_), boundary_(other.boundary_), 
				current_vertex_hedge_id_(other.current_vertex_hedge_id_), current_hedge_(other.current_hedge_)  {}
			
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S>& CirculatorTS<KERNEL,T,S>::operator=( const CirculatorTS<KERNEL,T,S>& rhs ) {
				current_hedge_ = rhs.current_hedge_;
				current_vertex_hedge_id_ = rhs.current_vertex_hedge_id_;
				boundary_ = rhs.boundary_;
				hedges_ = rhs.hedges_;
				return *this;
			}
			
			template <typename KERNEL, typename T, typename S>
			typename CirculatorTS<KERNEL,T,S>::reference CirculatorTS<KERNEL,T,S>::operator*( void ) {
				return getTS<KERNEL>(current_hedge_, dummyT, dummyS);
			}

			template <typename KERNEL, typename T, typename S>
			typename CirculatorTS<KERNEL,T,S>::pointer CirculatorTS<KERNEL,T,S>::operator->( void ) {
				return &getTS<KERNEL>(current_hedge_, dummyT, dummyS);
			}
			
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S>& CirculatorTS<KERNEL,T,S>::operator++( void ) {
				// Prima passo al successivo, poi restituisco il valore (incrementato)
				nextTS<KERNEL>(current_hedge_, hedges_, current_vertex_hedge_id_, boundary_, dummyT, dummyS);
				return (*this);
			}
			
			template <typename KERNEL, typename T, typename S>
			CirculatorTS<KERNEL,T,S> CirculatorTS<KERNEL,T,S>::operator++( int ) {
				// Prima passo al successivo, poi restituisco il valore (incrementato)
				CirculatorTS<KERNEL,T,S> temp(*this);
				++(*this);
				return temp;
			}
		}
	}
}

#endif
