#ifndef _MESH_MESHREADER_HPP
#define _MESH_MESHREADER_HPP

#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>

namespace ConservationLaw2D {
	namespace Mesh {
		/*! \namespace IO
		\brief Namespace per l'Input-Output per la mesh
		*/
		namespace IO {
			using namespace std;
			/*! \brief Legge una mesh da file*/
			template <typename MESH>
			void MeshReader( MESH& mesh, const string& filename ) {
				typedef typename MESH::real_t real_t;
				typedef typename MESH::Vertex_ptr vertex_ptr;
				typedef typename MESH::Polygon_ptr polygon_ptr;
				
				ifstream f(filename);
				if (f.fail()) {
					// Errore nell'apertura del file
					cout << "Error: meshfile '" << filename << "' not found!" << endl;
					exit(-1);
				}
				string currLine;
				size_t nV, nP, nE;
				// Leggo il numero di vertici e di lati
				do {
					getline(f, currLine);
				} while (currLine.find("# DATA") == string::npos);
				f >> nV >> nP >> nE;
				// VERTICI
				do {
					getline(f, currLine);
				} while (currLine.find("# POINTS") == string::npos);
				// Inserisco i vertici
				real_t x, y;
				vector<vertex_ptr> vhandle;
				vhandle.reserve(nV);
				for (size_t i = 0; i < nV; ++i) {
					f >> x >> y;
					vhandle[i] = mesh.addVertex(x, y);
				}
				// POLIGONI
				do {
					getline(f, currLine);
				} while (currLine.find("# ELEMENTS") == std::string::npos);
				// Inserisco i poligoni
				size_t nsides, id, color;
				polygon_ptr tmpp;
				vector<vertex_ptr> poly_vertex;
				poly_vertex.reserve(4);
				for (size_t i = 0; i < nP; ++i) {
					f >> nsides;
					for (size_t j = 0; j < nsides; ++j) {
						f >> id;
						assert(id < nV);
						poly_vertex.push_back( vhandle[id] );
					}
					f >> color;
					tmpp = mesh.addPolygon( poly_vertex );
					tmpp->setColor(color);
					poly_vertex.clear();
				}
				// LATI DI BORDO
				do {
					getline(f, currLine);
				} while (currLine.find("# EDGES") == std::string::npos);
				size_t v1, v2;
				bool found;
				typedef typename MESH::Vertex::HEdgeCirculator HECirc;
				for (size_t i = 0; i < nE; ++i) {
					f >> v1 >> v2 >> color;
					//cout << "Cerco lato " << *vhandle[v1] << "-" << *vhandle[v2] << " " << endl;
					// Cerco il lato di vertici v1 e v2
					HECirc vhec = vhandle[v1]->beginE();
					do {
						++vhec;
						//cout << " - " << *vhec << endl;
						found = ( &(vhec->vertexE()) == vhandle[v2] ) || ( &(vhec->vertexS()) == vhandle[v2] );
					} while( !found && vhec != vhandle[v1]->beginE() );
					if (found) {
						//cout << "  ==> Trovato!" << endl;
						vhec->setColor(color);
					} else {
						cout << "Lato " << *vhandle[v1] << "-" << *vhandle[v2] << " " << endl;
						cout << "  !!! Nessuna corrispondenza !!!" << endl;
						cout << " - " << v1 << endl;
						vhec = vhandle[v1]->beginE();
						++vhec;
						cout << " - " << *vhec << endl;
					}
				}
			} //MeshReader
		} // IO
	}
}

#endif
