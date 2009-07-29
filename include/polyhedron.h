#ifndef __POLYHEDRON_H__
#define __POLYHEDRON_H__

//
//  Polyhedron.h
//
//  Created by Robert Snapp on 2007-02-11.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include <algorithm>
#include <list>
#include <set>
#include <vector>
#include <iterator>
#include "ntuple.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <valarray>
#include "env.h"


#ifdef DEBUG
#define debug(str) std::err << str << std::endl
#else
#define debug(str)
#endif

#define printIndex(i) (i)

const size_t D = 3;

typedef double Point3d[D];
typedef float  Point3f[D];
typedef long   Point3i[D];
typedef size_t index_t;  // index_t will be used to refer to edge, vertex, and face indices.

// A polyhedron is an ordered set of three-dimensional vertices, edges, and 
// faces. The typename T designates the type of each vertex component (e.g.,
// char, short, int, float, or double.)
template<typename T>
class Polyhedron {
public:
  // typedef double normalType;
  typedef ntuple<T,D>  FaceVertex;
  typedef ntuple<T,D>  FaceNormal;
  typedef std::set<index_t, std::greater<index_t> >  IntSet;

  static const int SEARCH_FAILED = -1;
	
   Polyhedron() {} // default constructor
	
  // copy constructor
  Polyhedron(Polyhedron<T> const &p) {
    d_edgeTable   = p.d_edgeTable;
    d_vertexTable = p.d_vertexTable;
    d_faceTable   = p.d_faceTable;
	d_vMin = p.d_vMin;
	d_vMax = p.d_vMax;
  }
	
  // virtual destructor
  virtual ~Polyhedron() {}

  // addFace adds a new face, defined by the given vertex list, to the 
  // current polyhedral model. The index of the created face is returned.
  bool addFace(std::list<FaceVertex> const &vList, index_t* faceIndex);
	
  bool addFace(FaceVertex const &a, FaceVertex const &b, FaceVertex const &c, index_t* faceIndex) {
    std::list<FaceVertex> face;
    face.push_back(a);
    face.push_back(b);
    face.push_back(c);
    return addFace(face, faceIndex);	
  }
	
  // addFace(index_t, index_t, index_t) adds a triangular face with the indicated
  // indexed vertices. The index of the created face is returned.
  bool addFace(index_t v1Index, index_t v2Index, index_t v3Index, index_t* faceIndex);

  bool addFace(std::list<index_t> const &vIndexList, index_t* faceIndex);
	
  // addVertex adds the indicated vertex to the vertex table and associates it with the indicated edge,
  // or UNDEFINED if no edgeindex is provided. The function returns the vertexTable index of the new vertex.
  index_t addVertex(FaceVertex const &v, index_t edgeIndex = UNDEFINED);

  // marks the indicated face (along with its implicit edges, and vertices,
  // if they are unique) for deletion from the current model.
  // deleteFace returns true if the face was deleted without detecting any
  // errors. If false is returned, then an error was detected, and the 
  // calling routine should either throw an exception or abort.
  bool deleteFace(index_t faceIndex);
  bool deleteFace2(index_t faceIndex);
	
  // marks the indicated edge for deletion.
  bool deleteEdge(index_t edgeIndex);
  bool deleteEdge2(index_t edgeIndex);
	
  // marks the indicated vertex for deletion.
  bool deleteVertex(index_t vertexIndex);
  bool deleteVertex2(index_t vertexIndex);
	
  // Extends the face that contains the adjEdgeIndex by adding edges between 
  // the vertex with index v1Index and the endpoints of the edge adjEdgeIndex.
  // If successful, the boolean value true is returned. Otherwise, false
  // is returned.
  bool extendFace(index_t v1Index, index_t adjEdgeIndex);
	
  // getCollocalEdge searches for an edge different from edgeIndex that shares its vertex.
  // If found, the index of the discovered collocal edge is returned via the pointer in
  // the second argument, and the function returns the value true. If no edge can be
  // found, then false is returned. Note that three search strategies are employed in
  // succession.
  bool getCollocalEdge(index_t edgeIndex, index_t *collocalEdge) const;

  ntuple<T,D> getMax() {
	return d_vMax;
  }

  ntuple<T,D> getMin() {
	return d_vMin;
  } 
	
  // removes all faces, halfedges, and vertices that are marked for deletion.
  // purgeFaces returns the number of faces that were purged.
  bool purgeFaces();
  bool purgeEdges();
  bool purgeVertices();
  bool purgeTables();
		
	
  // Returns the number of faces in the current polytope.
  index_t faces() const    {return d_faceTable.size();  }
	
  // Returns the number of edges in the current polytope.
  index_t edges() const    {return d_edgeTable.size();  }
	
  // Returns the number of vertices in the current polytope.
  index_t vertices() const {return d_vertexTable.size();}
	
  // isConvex() returns true of the polyhedron is convex.
  bool isConvex() const;
	
  // isColinear(vIndex, eIndex) returns true if the specified vertex and edge are
  // colinear.
  bool isColinear(index_t vIndex, index_t eIndex) {
    ntuple<T,D> v1 = d_vertexTable[vIndex].d_v;
    ntuple<T,D> v2 = d_vertexTable[d_edgeTable[eIndex]._v];
    ntuple<T,D> v3 = d_vertexTable[d_edgeTable[d_edgeTable[eIndex]._next]];
    T tp = tripleProduct(v1, v2, v3); 

    return tp != 0;
  }
  // Returns the next available index for the vertex table.
  index_t nextVertexIndex() {return vertices();}
	
  // Returns the next available index for the edge table.
  index_t nextEdgeIndex() {return edges();}
	
  // Returns the next available index for the edge table.
  index_t lastEdgeIndex() {return edges() - 1;}
	
  // Returns the next available index for the face table.
  index_t nextFaceIndex() {return faces();}
	
  // Returns the number of edges that lack a twinEdge.
  int boundaryEdges() const {
    int count = 0;
    for(index_t i = 0; i < d_edgeTable.size(); ++i) {
      if (d_edgeTable[i].d_twinEdge == UNDEFINED) 
        ++count;
    }
    return count;
  }
	
  FaceNormal &getNormalVector(index_t faceIndex) {
    assert(0 <= faceIndex && faceIndex < faces());
    return d_faceTable[faceIndex].d_normalVector;
  }
	
  FaceVertex &getFaceVertex(index_t faceIndex) {
    assert(0 <= faceIndex && faceIndex < faces());
    return d_vertexTable[d_edgeTable[d_faceTable[faceIndex].d_outer].d_vertex].d_v;
  }

#ifdef COMMENT
  // A vertex v is said to be exterior to face f, if v lies on the
  // positive side of the plane that coincides with f. Vertex v, is said to be
  // coplanar with respect to f, if v lies in the plane that coincides with f,
  // and is said to be interior to face f, if v lies on the negative side
  // of the plane that coincides with face f. The function visibilityTest(...) 
  // returns +1 if the vertex v is exterior, 0 if it is coplanar, and -1 if v
  // is interior to f. For queries of floating point data, the third argument
  // epsilon can be used to define a slab of with epsilon that is centered
  // on face f.
  int visibilityTest(index_t faceIndex, FaceVertex& v, T epsilon=0) {
  }
#endif

  // Returns the number of edges that have a twinEdge.
  int interiorEdges() const {
    return edges() - boundaryEdges();
  }
	
  // Checks the three tables (d_faceTable, d_edgeTable, and d_vertexTable) for the presence of
  // inactive elements. If the first argument is set to true, then a detailed error report
  // is printed on std::cerr. If no inactive elements are detected, then checkTables returns
  // true, otherwise false.
  bool checkTables(bool verbose) const;

  // Function countInactiveFaces returns the number of faces in the d_faceTable that are not
  // active, i.e. have d_isPresent set to false.
  index_t countInactiveFaces() const;

  // Function countInactiveEdges returns the number of edges in the d_edgeTable that are not
  // active, i.e. have d_isPresent set to false.
  index_t countInactiveEdges() const;

  // Function countInactiveVertices returns the number of vertices in the d_vertexTable that are not
  // active, i.e. have d_isPresent set to false.
  index_t countInactiveVertices() const;

  // Function eulerTest tests the topological consistency of the polyhedron defined
  // by the three tables: d_faceTable, d_edgeTable, and d_vertexTable. If the argument
  // verbose is true, then a detailed error report is printed on std::cerr. If no
  // topological errors are detected, then eulerTest returns true, otherwise false.
  bool eulerTest(bool verbose = false) const; 
	
  // Prints a description of the current polyhedron on the indicated stream.
  void describe(std::ostream &os) const;
	
  // Computes are returns the area of the face with index i.
  double faceArea(index_t i) const;
 
  // Computes and returns the surface are of the current polyhedron.
  double surfaceArea() const ;

  // Computes and returns the volume of the current polyhedron.
  double volume() const;

  // Display the polyhedron using OpenGL
  void display() const;

  // friend std::ostream& operator<<(std::ostream &os, Polyhedron<T>::HalfEdge const &e);
protected:
  class HalfEdge;  // Nested class for representing each HalfEdge
  class Vertex;    // Nested class for representing each Vertex
  class Face;      // Nested class for representing each Face
	
  // Serves as a null pointer.
  static const index_t UNDEFINED = static_cast<index_t>(-1);  // Set UNDEFINED to largest index_t.
  static const int poly_level = 1; // debug level
	
#if (DEBUG 	> 0)	
  static const bool debug = (DEBUG >= poly_level);
#else
  static const bool debug = false;
#endif
	
  // Compute the vertex extrema, that is the values of d_vMin and d_vMax, using the current values
  // stored in d_vertexTable. Normally, these values are updated by addVertex, and purgeVertices.
  // However, in the event that the vertex table has been modified by other means, this function
  // should be invoked.

  void computeVertexExtrema(void);

  index_t getUndefined() const {return UNDEFINED;}
	
  // Returns true if a halfedge already exists that extends from vertex p 
  // to vertex n, within precision epsilon. If true, then index returns the
  // location of this edge in the HalfEdge table.
  bool findEdge(FaceVertex const &p, FaceVertex const &n, index_t *index, 
                double epsilon = 0) const;
	
  // findEdge(index_t, index_t, index_t) returns true if there is an edge in the edge
  // table that extends from the initial vertex to the final vertex. In
  // this case, the index of the discovered edge is returned via the 
  // third argument. If no edge is found, false is returned.
  bool findEdge(index_t initialVertexIndex, index_t finalVertexIndex, index_t *index) const;
	
  // Returns true if a vertex v already exists in the vertex table, within 
  // precision epsilon. If true, then index returns the location of this 
  // vertex in the vertex table.
  bool findVertex(FaceVertex const &v, index_t* index, double epsilon = 0) const;
	
  // Returns the index of the current Face for the indicated edge.
  index_t getCurrentFace(index_t edge) const {
    return (isDefined(edge) ? d_edgeTable[edge].d_face : getUndefined());
  }
	
  // Returns the index of the current Vertex for the indicated edge.
  index_t getCurrentVertex(index_t edge) const {
    return (isDefined(edge) ? d_edgeTable[edge].d_vertex : getUndefined());
  }
	
  // Returns the halfedge index that follows the value edge.
  index_t getNextEdge(index_t edge) const {
    return (isDefined(edge) ? d_edgeTable[edge].d_nextEdge : getUndefined());
  }
	
	
  index_t getClockwiseEdge(index_t edge) const {
    return getNextEdge(getTwinEdge(edge));
  }
	
  index_t getCounterClockwiseEdge(index_t edge) const {
    return getTwinEdge(getPrevEdge(edge));
  }
	
  // Returns the index of the next Vertex for the indicated edge.
  index_t getNextVertex(index_t edge) const {
    return getCurrentVertex(getNextEdge(edge));
  }
	
  // Returns the previous half edge
  index_t getPrevEdge(index_t edge) const {
    return (isDefined(edge) ? d_edgeTable[edge].d_prevEdge : getUndefined());
  }
		
  // Returns the halfedge index that precedes the value edge.
  index_t getPrevVertex(index_t edge) const {
    return getCurrentVertex(getPrevEdge(edge));
  }
	
  // Returns the halfedge index that is opposite to the value edge.
  index_t getTwinEdge(index_t edge) const {
    return (isDefined(edge) ? d_edgeTable[edge].d_twinEdge : getUndefined()); 
  }
	
  // Returns true if the halfedge with index value edge has a twin.
  bool hasTwin(index_t edge) const {
    return (getTwinEdge(edge) != UNDEFINED);
  }
	
  bool isDefined(index_t index) const {
    return (index != UNDEFINED);
  }

  // Function isolatedVertexCheckAndDelete deletes the edge with value edgeToDelete after
  // checking to see if its vertex is used by any other edge in the edge circuit that 
  // contains the edge with index successor. If no other edge has the same vertex, then
  // the vertex is also deleted. Otherwise, the edge field of the vertex entry is updated
  // if necessary, to ensure that there exists a reference to an active edge. Also, the
  // outer edge field of the face entry is updated, if necessary. Note that edgeToDelete
  // and successor should belong (formerly) to the same face.
  bool isolatedVertexCheckAndDelete(index_t edgeToDelete, index_t successor);

  // Returns true if and only if the argument is the index for an edge that is present,
  // that is not deleted, in the current polyhedron. Note that if edge is UNDEFINED, then
  // false is returned.
  bool isPresentEdge(index_t index) const {
    return (0 <= index && index < edges() && d_edgeTable[index].d_isPresent);
  }
	
  index_t edgeOfVertex(index_t vertexIndex) {
    return (isDefined(vertexIndex) ? d_vertexTable[vertexIndex].d_edge : getUndefined());
  }
	
  bool updateEdgesOfVertex(index_t vIndex, index_t vIndexOld);
	
  bool vertexIsIsolated(index_t vertexIndex) {
    return edgeOfVertex(vertexIndex) == UNDEFINED; 
  }
	
  // Compute the normal vector to the indicated face.
  void computeNormalVector(index_t faceIndex);
	
  // Compute the normal vectors of every face.
  void computeNormalVectors();
	
  // Merges the two indicated faces if they are adjacent and coplanar.
  bool mergeFaces(index_t fold, index_t fnew);
	
  // Clear the polytope.
  void clear() {
    d_edgeTable.clear();
    d_vertexTable.clear();
    d_faceTable.clear();
  }
	
  // ================
  // = Data members =
  // ================
	
  std::vector<HalfEdge> d_edgeTable;
  std::vector<Vertex>   d_vertexTable;
  std::vector<Face>	    d_faceTable;

  IntSet    d_deadFaces;
  IntSet    d_deadEdges;
  IntSet    d_deadVertices;

  ntuple<T,D> d_vMax;	// maintains a three dimensional bounding box that contains
  ntuple<T,D> d_vMin;   // the polyhedral model
};

template<typename T>
class Polyhedron<T>::HalfEdge {
public:
  // Constructors
  explicit HalfEdge() : d_vertex(UNDEFINED), d_prevEdge(UNDEFINED), 
                        d_nextEdge(UNDEFINED), d_twinEdge(UNDEFINED), d_face(UNDEFINED), 
	                    d_isPresent(true) {}
	
  HalfEdge(index_t pVertex, index_t lpEdge, index_t lnEdge, index_t twinEdge, index_t lFace) : 
    d_vertex(pVertex), d_prevEdge(lpEdge), d_nextEdge(lnEdge), 
    d_twinEdge(twinEdge), d_face(lFace), d_isPresent(true) {}
	
  // Copy Constructor
  HalfEdge(const HalfEdge& he) : d_vertex(he.d_vertex), 
                                 d_prevEdge(he.d_prevEdge), d_nextEdge(he.d_nextEdge), 
                                 d_twinEdge(he.d_twinEdge),
                                 d_face(he.d_face), d_isPresent(he.d_isPresent) {}
	
  HalfEdge& operator=(const HalfEdge &he) {
    d_vertex   = he.d_vertex;
    d_prevEdge = he.d_prevEdge;
    d_nextEdge = he.d_nextEdge;
    d_twinEdge = he.d_twinEdge;
    d_face     = he.d_face;
    d_isPresent  = he.d_isPresent;
    return *this;
  }
	
  void setVertex(index_t pVertex) {
    d_vertex = pVertex;
  }
	
  void setPrevEdge(index_t lpEdge) {
    d_prevEdge = lpEdge;
  }
	
  void setNextEdge(index_t lnEdge) {
    d_nextEdge = lnEdge;
  }
	
  void setTwinEdge(index_t twinEdge) {
    d_twinEdge = twinEdge;
  }
	
  void setFace(index_t lFace) {
    d_face = lFace;
  }
	
  index_t getVertex()   const {return d_vertex; }
  index_t getPrevEdge() const {return d_prevEdge; }
  index_t getNextEdge() const {return d_nextEdge; }
  index_t getTwinEdge() const {return d_twinEdge;}
  index_t getFace()     const {return d_face; }
	
  bool operator==(HalfEdge& he) {
    return d_vertex   == he.d_vertex 
      && d_prevEdge == he.d_prevEdge 
      && d_nextEdge == he.d_nextEdge 
      && d_twinEdge == he.d_twinEdge 
      && d_face     == he.d_face
      && d_isPresent  == he.d_isPresent;
  }

#ifdef COMMENT  // TODO: Repair output functions
  std::ostream& operator<< <>(std::ostream &os, HalfEdge &e) {
    std::ostringstream s;
    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision);

    s << "HalfEdge(" 
      << e.d_face << ", " 
      << e.d_prevEdge << ", " 
      << e.d_nextEdge << ", " 
      << e.d_twinEdge << ", " 
      << e.d_vertex << ","
      << e.d_isPresent << ")";
    os << s.str();
    return os;
  }
#endif

  index_t d_vertex; 	
  index_t d_prevEdge;  	
  index_t d_nextEdge;  
  index_t d_twinEdge; 	
  index_t d_face;	
  bool d_isPresent;  	// true if the Face exists; false if deleted;
};

template<typename T>
class Polyhedron<T>::Face {
public:

	
  Face() : d_outer(UNDEFINED), d_isPresent(true) {
    d_normalVector.clear();
  }
	
  explicit Face(index_t outerEdge) : d_outer(outerEdge), d_isPresent(true) {
    d_normalVector.clear();
  }
	
  void setOuterEdge(index_t e) {d_outer = e;}
  
  index_t 				d_outer;
  std::vector<index_t> 	d_inner;	
  FaceNormal            d_normalVector;
  bool 				    d_isPresent;  	// true if the Face exists; false if deleted;
};


template<typename T>
class Polyhedron<T>::Vertex {
public:

  Vertex() : d_edge(UNDEFINED), d_isPresent(true) {
    d_v.clear();
  }
	
  Vertex(Point3d &v, index_t edge) : d_edge(edge), d_isPresent(true) {
    d_v.resize(D);
    for (index_t i = 0; i < D; i++) {
      d_v[i] = static_cast<T>(v[i]);
    }
  }
	
  Vertex(FaceVertex const &v, index_t const edge) : d_v(v) , d_edge(edge), d_isPresent(true) {}
	
	
  // copy constructor
  Vertex(Vertex const &v) {
    d_v    = v.d_v;
    d_edge = v.d_edge;
    d_isPresent = v.d_isPresent;
  };
	
  // withinEpsilon should return true if and only if the vertex is located 
  // within a distance epsilon of the position w, using an Lp metric indexed
  // by p. (The default metric is Euclidean.)	
  bool isWithinEpsilon(	FaceVertex const &w, 
                        double const epsilon, 
                        double const p=2.0) const {
    return d_v.isWithinEpsilon(w, epsilon, p);
  }
	
  bool isWithinEpsilon(	Polyhedron<T>::Vertex const &w, 
                        double const epsilon, 
                        double const p=2.0) const {
    return d_v.isWithinEpsilon(w.d_v, epsilon, p);
  }
				
  FaceVertex 	d_v;		  // An n-tuple of type T.
  index_t 			d_edge;    // HalfEdge* d_edge;
  bool			d_isPresent;
};


// =====================
// = Private Functions =
// =====================

// TODO: improve the efficiency of the following function by searching through the
// vertex table for p, and then use the topology of the edgeTable to find n.
template<typename T>
inline bool Polyhedron<T>::findEdge(FaceVertex const &p, 
									FaceVertex const &n, 
									index_t *index, 
									double epsilon) const {
  for(index_t i = 0; i < edges(); i++) {
    if (d_vertexTable[getCurrentVertex(i)].isWithinEpsilon(p, epsilon) && 
        d_edgeTable[i].d_isPresent) {
      index_t theNextVertexIndex = getNextVertex(i);
      if (isDefined(theNextVertexIndex)) {
        if (d_vertexTable[theNextVertexIndex].isWithinEpsilon(n, epsilon)) {
          *index = i;
          return true;
        }
      }
    }
  }
  return false;
}

// There is a problem in the following.

template<typename T>
inline bool Polyhedron<T>::findEdge(	index_t initialVertexIndex, 
										index_t finalVertexIndex, 
										index_t *index) const {
	// ensure that the first two arguments point to existing vertices.
	assert(0 <= initialVertexIndex && initialVertexIndex < vertices());
	assert(0 <= finalVertexIndex   && finalVertexIndex   < vertices());
	
	// retain the index of the HalfEdge associated with the initial vertex.
	index_t initialEdgeIndex = d_vertexTable[initialVertexIndex].d_edge;
	
	if (! isPresentEdge(initialEdgeIndex)) {
		// std::cerr << "Polyhedron<T>::findEdge: initial edge for vertex " << initialVertexIndex
		// 			<< " is not present in the edge table." << std::endl;
		return false;
	}
	
	// As we walk around the vertex (clockwise), edge indicates the current
	// edge that points away from the initial vertex.
	//	index_t nEdges = edges();
	index_t edge = initialEdgeIndex;
	do {
		if (getNextVertex(edge) == finalVertexIndex) {  // Eureka!
			*index = edge; // return the discovered edge index.
			return true;
		}
		
		// Search the edges that eminanate from the initial vertex in a clockwise order.
		edge = getClockwiseEdge(edge);
	} while(edge != initialEdgeIndex && isPresentEdge(edge));
	
	if (edge != initialEdgeIndex) {	
	  // Starting again from the initialEdge for the initial vertex, search the edges that
	  // eminate from the initial vertex in a counter-clockwise order.
		
	  // since the initial edge was checked previously, skip forward.
	  initialEdgeIndex = getCounterClockwiseEdge(initialEdgeIndex);  
		
	  if (isPresentEdge(initialEdgeIndex)) {
			index_t edge = initialEdgeIndex;
			do {
				if (getNextVertex(edge) == finalVertexIndex) {  // Eureka!
					*index = edge; // return the discovered edge index.
					return true;
				}

				edge = getCounterClockwiseEdge(edge); // search counter-clockwise about initial vertex.
			
				// Continue walking until either the twin is undefined, or we
				// return to the initalEdgeIndex.
			} while(edge != initialEdgeIndex && isPresentEdge(edge));	
		}
	}
	
	// Now search around the final vertex
	initialEdgeIndex = getPrevEdge(d_vertexTable[finalVertexIndex].d_edge);
	if (! isPresentEdge(initialEdgeIndex)) {
		initialEdgeIndex = getTwinEdge(d_vertexTable[finalVertexIndex].d_edge);
	}
	
	if (isPresentEdge(initialEdgeIndex)) {
		edge = initialEdgeIndex;
		do {
			if (getCurrentVertex(edge) == initialVertexIndex) {  // Eureka!
				*index = edge; // return the discovered edge index.
				return true;
			}
		
			edge = getTwinEdge(getNextEdge(edge));	
			// Continue walking until either the twin is undefined, or we
			// return to the initalEdgeIndex.
		} while(edge != initialEdgeIndex && isPresentEdge(edge));
	
		if (edge != initialEdgeIndex) {
			// Search for a matching edge in the counter-clockwise direction about the final vertex.
			index_t edge = getPrevEdge(getTwinEdge(initialEdgeIndex)); // counter-clockwise step
			if (isPresentEdge(edge)) {
				do {
					if (getCurrentVertex(edge) == initialVertexIndex) {  // Eureka!
						*index = edge; // return the discovered edge index.
						return true;
					}

					edge = getPrevEdge(getTwinEdge(edge)); // counter-clockwise step
					// Continue walking until either the twin is undefined, or we
					// return to the initalEdgeIndex.
				} while(edge != initialEdgeIndex && isPresentEdge(edge));	
			}
		}
	}	
	return false;  // Twin edge was not found.
}									

template<typename T>
inline bool Polyhedron<T>::findVertex(	FaceVertex const &v, 
										index_t *index, 
										double epsilon) const {	
	for(index_t i = 0; i < vertices(); i++) {
		if (d_vertexTable[i].isWithinEpsilon(v, epsilon)) {
			*index = i;
			return true;
		}
	}
	return false;
}


// ====================
// = Public Functions =
// ====================

// addVertex adds a new vertex record to the vertex table. The second argument should
// be assigned to the index of the half edge rooted at the new vertex, or UNDEFINED.
// The index of the new vertex record is returned.  addVertex also updates the bounding
// box coordinates, d_vMin and d_vMax.
template<typename T>
inline index_t Polyhedron<T>::addVertex(FaceVertex const &v, index_t edgeIndex) {
  Vertex newVertex(v, edgeIndex);
  if (vertices() == 0) {
	d_vMax = newVertex.d_v;
	d_vMin = newVertex.d_v;
  } else {
	for(size_t i = 0; i < D; i++) {
	  T x = newVertex.d_v[i];
	  d_vMax[i] = std::max<T>(d_vMax[i], x);
	  d_vMin[i] = std::min<T>(d_vMin[i], x);
	}
  }
  d_vertexTable.push_back(newVertex);
  return d_vertexTable.size() - 1;
}

// addFace of a list of ntuples provides a convenient interface for defining a new face record
// in the face table that coincides with the polygon defined by this coordinate list. The index
// of the new face is returned through the second argument. Adding a face is somewhat complicated
// as new records probably need to be added to both the vertex and edge tables as well, and existing
// fields may need updating, including twin edge indicies. Most of the heavy lifting is performed
// in the function addFace that accepts a list of vertex indices as its initial argument. Thus,
// the following routine simply ensures that a record for each vertex in the vertex list v exists
// in the vertex table, creating new vertex records as needed, and then invokes the other version
// of addFace using a constructed list of vertex indices. If all succeeds, the value true is returned;
// if not, the function returns false.
template<typename T>
inline bool Polyhedron<T>::addFace(typename std::list< FaceVertex > const &v, index_t* faceIndex) {
  typename std::list< FaceVertex >::const_iterator preVertex;
  std::list<index_t> vIndicies;
  
   for(preVertex = v.begin(); preVertex != v.end(); ++preVertex) {
     index_t vertexIndex;
     if (! findVertex(*preVertex, &vertexIndex)) {
       // append *preVertex to the vertex table, with an UNDEFINED edge.
	   vertexIndex = addVertex(*preVertex);
     }
     vIndicies.push_back(vertexIndex);
   }

   return Polyhedron<T>::addFace(vIndicies, faceIndex);
}

// The following implementation of addFace will construct a new face in the polyhedron that is delimitted
// by the vertices that correspond to the list of vertex indices in the first arguement. The index of the
// new face in the face table is stored in *faceIndex. The member function returns true if the operation
// was successful.

template<typename T>
bool Polyhedron<T>::addFace(std::list<index_t> const &vIndexList, index_t* faceIndex) {
  typename std::list< index_t >::const_iterator preVertex;
  index_t size = vIndexList.size(); // the number of vertices and edges in the new face.
  index_t edgeIndexBase = edges();  // the index of the next new edge.
  //index_t edgeIndex = edgeIndexBase;
  *faceIndex = faces();             // since this face will be appended to the face table,
                                    // the index of the new face is the current face table size.
  index_t count = 0;

  Face newFace(edgeIndexBase);
  d_faceTable.push_back(newFace);    // Optimism!
  
  // Cycle through the vertices of the new face. Note that edgeIndex denotes the current edge.
  for(preVertex = vIndexList.begin(); preVertex != vIndexList.end(); ++preVertex, //++edgeIndex,
     ++count) {
    assert(*preVertex >= 0 && *preVertex < vertices());  // Ensure that the vertex indices are valid.

	index_t edgeIndex     = edgeIndexBase + count;

	// It is a mistake to insert a reference to an edge with an index greater than
	// or equal to edgeIndexBase + count, until that particular HalfEdge has been
	// added to the edge table, as function findEdge called below might try to
	// access a nonexistent entry in the edge table. Thus, for the time being, we set 
	// these references to UNDEFINED.

	//	index_t nextEdgeIndex = (count == size - 1 ? edgeIndexBase : edgeIndex + 1);
	index_t nextEdgeIndex = (count == size - 1 ? edgeIndexBase : UNDEFINED);
	index_t prevEdgeIndex = (count == 0 ? UNDEFINED : edgeIndex - 1);

    // Find the twin edge, if it exists. First step: find the next vertex
    // in the list of points.
    typename std::list< index_t >::const_iterator nextVertex = preVertex;
    if (distance(preVertex, vIndexList.end()) > 1) {
      ++nextVertex; 
    } else {
      // for the last edge, the next vertex is the initial one.
      nextVertex = vIndexList.begin();
    }

     // The twin edge should point from the next vertex to the current one.
    index_t twinEdgeIndex;
    if (findEdge(*nextVertex, *preVertex, &twinEdgeIndex)) {
      if (isDefined(d_edgeTable[twinEdgeIndex].d_twinEdge)) {
        std::cerr << 
          "Polyhedron<T>::addFace: ERROR: current face may conincide with a previous face."
                  << std::endl;
		// This return is awkward, as we have not yet removed any half edges or vertices from the
		// current (illegal) face. Three fixes come to mind. Either allow twin edges to degenerate,
		// as one might have if the (noncovex) polyhedron folds back on itself; or add half edges
		// initially with the present flag set to false; then once all is well, set them all to
		// true; or go through the face again up to the point of failure, removing all new halfedges
		// from edgeIndexBase and up.
        return false;
      }
      d_edgeTable[twinEdgeIndex].d_twinEdge = edgeIndex; // This should be okay, as it's after findEdge.
    } else {
      twinEdgeIndex = UNDEFINED;
    }

     HalfEdge newEdge(*preVertex, prevEdgeIndex, nextEdgeIndex, 
                     twinEdgeIndex, *faceIndex);
	 d_edgeTable.push_back(newEdge);

	 // Now it is safe to refere to edgeIndex in the edgeTable.
	 
	 if (count > 0) {
	   assert(! isDefined(d_edgeTable[edgeIndex - 1].d_nextEdge));
	   d_edgeTable[edgeIndex - 1].d_nextEdge = edgeIndex;
	 }

	 // If the current vertex is isolated, then connect it to the current edge.
	 if (! isDefined(d_vertexTable[*preVertex].d_edge)) { 
	   d_vertexTable[*preVertex].d_edge = edgeIndex;
	 }
     
	 // std::cerr << "Added edge " << edges() - 1 << " from vertex " << *preVertex 
     //           << " to vertex " << *nextVertex << std::endl;
  }
  
  assert(! isDefined(d_edgeTable[edgeIndexBase].d_prevEdge));
  d_edgeTable[edgeIndexBase].d_prevEdge = edgeIndexBase + size - 1;
    
  computeNormalVector(*faceIndex);
  return true;
}

template<typename T>
inline bool Polyhedron<T>::addFace(index_t v1Index, index_t v2Index, index_t v3Index, index_t* faceIndex) {
	assert(0 <= v1Index && v1Index < vertices());
	assert(0 <= v2Index && v2Index < vertices());
	assert(0 <= v3Index && v3Index < vertices());
	
	index_t twin21 = UNDEFINED;
	index_t twin13 = UNDEFINED;
	index_t twin32 = UNDEFINED;
		
	index_t nextEdgeIndex = edges();
	*faceIndex = faces();
	
	if (findEdge(v2Index, v1Index, &twin21)) {
		index_t oldTwin = d_edgeTable[twin21].d_twinEdge;
		if (isDefined(oldTwin) && d_edgeTable[oldTwin].d_isPresent) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin21
						<< " connecting vertex " << v2Index << " to " << v1Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
		//d_edgeTable[twin21].d_twinEdge = nextEdgeIndex;
	} else {
		twin21 = UNDEFINED;
	}
	
	if (findEdge(v3Index, v2Index, &twin32)) {
		index_t oldTwin = d_edgeTable[twin32].d_twinEdge;
		if (isDefined(oldTwin) && d_edgeTable[oldTwin].d_isPresent) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin32
						<< " connecting vertex " << v3Index << " to " << v2Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
		//d_edgeTable[twin32].d_twinEdge = nextEdgeIndex + 1;
	} else {
		twin32 = UNDEFINED;
	}
	
	if (findEdge(v1Index, v3Index, &twin13)) {
		index_t oldTwin = d_edgeTable[twin13].d_twinEdge;
		if (isDefined(oldTwin) && d_edgeTable[oldTwin].d_isPresent) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin13
						<< " connecting vertex " << v1Index << " to " << v3Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
		// d_edgeTable[twin13].d_twinEdge = nextEdgeIndex + 2;
	} else {
		twin13 = UNDEFINED;
	}
	
	if (isDefined(twin21)) d_edgeTable[twin21].d_twinEdge = nextEdgeIndex;
	if (isDefined(twin32)) d_edgeTable[twin32].d_twinEdge = nextEdgeIndex + 1;
	if (isDefined(twin13)) d_edgeTable[twin13].d_twinEdge = nextEdgeIndex + 2;
	
	if (vertexIsIsolated(v1Index)) {
		d_vertexTable[v1Index].d_edge = nextEdgeIndex;
	}
	
	if (vertexIsIsolated(v2Index)) {
		d_vertexTable[v2Index].d_edge = nextEdgeIndex + 1;
	}
	
	if (vertexIsIsolated(v3Index)) {
		d_vertexTable[v3Index].d_edge = nextEdgeIndex + 2;
	}
	
	HalfEdge edge12(v1Index, 
					nextEdgeIndex + 2, nextEdgeIndex + 1, twin21, 
					*faceIndex);
	HalfEdge edge23(v2Index, 
					nextEdgeIndex,     nextEdgeIndex + 2, twin32, 
					*faceIndex);
	HalfEdge edge31(v3Index, 
					nextEdgeIndex + 1, nextEdgeIndex,     twin13, 
					*faceIndex);
					
	d_edgeTable.push_back(edge12);
	d_edgeTable.push_back(edge23);
	d_edgeTable.push_back(edge31);
	
	Face newFace(nextEdgeIndex);
	d_faceTable.push_back(newFace);
	computeNormalVector(*faceIndex);
	return true;
}

// Function isolatedVertexCheckAndDelete is an auxilliary function used by function extendFace
// (see below).The first argument, edgeIndex, should identify a half-edge (E) that belongs to a
// closed polygonal face (F) contained in the internal tables. If the vertex (V) of this half-edge
// does not belong to any other half-edge in F, then both the V and E will be deleted by inserting
// them in the d_deadVertices and d_deadEdges sets respectively. If V belongs to another half-edge,
// which can occur if the edges of the face re-enter V (like a simple hour glass), then only the
// E is deleted. In this case, the edge that belongs to V will be changed if necessary, to an
// active edge index. 
//
// The second argument, goodEdge, should evaluate to an active edge index in the current face.
// This argument is provided as pointers to edgeToDelete from its antecedent or successor, may
// have been removed. Thus the nextEdge links from goodEdge and its successors should form a
// closed traversal of a face that may no longer contain edgeToDelete.
//
// If all edges of the face are defined, and form a closed loop, then the above is performed
// and the value true is returned. The value false is returned if the face is incomplete, in
// which case a calling function should either throw an exception or abort the process.

template<typename T>
bool Polyhedron<T>::isolatedVertexCheckAndDelete(index_t edgeToDelete, index_t goodEdge) {
  index_t vIndex = getCurrentVertex(edgeToDelete);
  index_t altEdge = UNDEFINED;
  index_t edge = goodEdge;

  // Traverse the face inorder to find another edge that references vIndex.
  do {
    if (! isDefined(edge)) {
      std::cerr << "Polyhedron<T>::isolatedVertexCheckAndDelete(" 
                << edgeToDelete << ", " << goodEdge << "): "
                << "Undefined edge encounted." << std::endl;
      return false;
    }

    if (getCurrentVertex(edge) == vIndex && edge != edgeToDelete) {
      altEdge = edge;
      if (d_vertexTable[vIndex].d_edge == edgeToDelete) {
        d_vertexTable[vIndex].d_edge = altEdge;
      }
      break;
    }
          
    edge = getNextEdge(edge);
  } while (edge != goodEdge);

  if (! isDefined(altEdge)) {   // No edge was found
    d_deadVertices.insert(vIndex);
    d_vertexTable[vIndex].d_isPresent = false;
  }

  index_t faceIndex = getCurrentFace(edgeToDelete);
  if (d_faceTable[faceIndex].d_outer == edgeToDelete) {
    d_faceTable[faceIndex].d_outer = goodEdge;
  }

  d_deadEdges.insert(edgeToDelete);
  d_edgeTable[edgeToDelete].d_isPresent = false;
  return true;
}


// Function extendFace adds a coplanar triangular region (R) to an existing face (F). Two
// arguments are provided that define the R and F: v1Index is the index of vertex V1 that 
// belongs to R, and adjEdgeIndex is the index of an edge, directed from vertex V3 to V2,
// that belongs to F (but does not contain V1). Thus R is defined as the triangle that 
// contains both adjEdgeIndex and V1, hence, triangle V1-V2-V3.
//
// Since F need not be convex (which may occur during the construction phase) four distinct
// cases are possible:
//
// 1. V1 does not belong to F (point extension);
// 2. Edge V2->V1 belongs to F, but edge V1->V3 does not (left notch fill);
// 3. Edge V1->V3 belongs to F, but edge V2->V1 does not (right notch fill);
// 4. Edges V2->V1 and V1->V3 belong to F, (hole fill).
//
//  Case 1 requires the creation of one new edge (v1->v2).
//  Case 2 requires the deletion of one edge (v2->v1) and vertex v2, provided it is not
// visited by a re-entrant edge. The adjacent edge (v3->v2) is displaced to form v3->v1.
//  Case 3, likewise, requires the deletion the adjacent Edge (v3->v2), and vertex v3, provided
// it is not visited by a re-entrant edge. The existing edge (v1->v3) is displaced to form v1->v2.
//  Case 4 requires the deletion of all three edges (v1->v3, v3->v2, and v2->v1).

template<typename T>
inline bool Polyhedron<T>::extendFace(index_t v1Index, index_t adjEdgeIndex) {
	assert(0 <= v1Index && v1Index < vertices());               // Validate arguments.
	assert(0 <= adjEdgeIndex && adjEdgeIndex < edges());
	
	index_t twin21;                                                 // v2->v1
	index_t twin13;                                                 // v1->v3
	index_t v2Index   = getNextVertex(adjEdgeIndex);                // v2
	index_t v3Index   = getCurrentVertex(adjEdgeIndex);             // v1
	index_t faceIndex = getCurrentFace(adjEdgeIndex);               // F
	

	if (findEdge(v2Index, v1Index, &twin21)) {                  // Does v2->v1 exist?
		index_t oldTwin = d_edgeTable[twin21].d_twinEdge;             // If so, it should be twinless.
		if (isDefined(oldTwin) && d_edgeTable[oldTwin].d_isPresent) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin21
						<< " connecting vertex " << v2Index << " to " << v1Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
	} else {
		twin21 = UNDEFINED;
	}
	
	if (findEdge(v1Index, v3Index, &twin13)) {                  // Does v1->v3 exist?
		index_t oldTwin = d_edgeTable[twin13].d_twinEdge;             // If so, it should be twinless.
		if (isDefined(oldTwin) && d_edgeTable[oldTwin].d_isPresent) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin13
						<< " connecting vertex " << v1Index << " to " << v3Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
	} else {
		twin13 = UNDEFINED;
	}
	
	if (getCurrentFace(twin21) == faceIndex) {                  // Case 2 or 4?
		index_t successor = getNextEdge(twin21);
    
        if (! isolatedVertexCheckAndDelete(twin21, successor)) {
          std::cerr << "Polyhedron<T>::extendFace(" << v1Index << ", "
                    << adjEdgeIndex << "): " 
                    << "Undefined edge encounted while traversing face "
                    << faceIndex << std::endl;
          return false;
        }

		if (getCurrentFace(twin13) == faceIndex) {              // Case 4?
			// The extension fills in a triangular hole.
			index_t antecedent = getPrevEdge(twin13);
			d_edgeTable[antecedent].d_nextEdge = successor;
			d_edgeTable[successor ].d_prevEdge = antecedent;
              
            deleteEdge(twin13);
            
            if (! isolatedVertexCheckAndDelete(adjEdgeIndex, successor)) {
              std::cerr << "Polyhedron<T>::extendFace(" << v1Index << ", "
                        << adjEdgeIndex << "): " 
                        << "Undefined edge encounted while traversing face "
                        << faceIndex << std::endl;
              return false;
            }

		} else {                                                // Case 2.
			// The extension fills a notch between vertices v3, v2, v1.
			d_edgeTable[adjEdgeIndex].d_nextEdge = successor;
			d_edgeTable[adjEdgeIndex].d_twinEdge = twin13;
			if (isDefined(twin13)) {
				d_edgeTable[twin13].d_twinEdge = adjEdgeIndex;
			}
			d_edgeTable[successor].d_prevEdge = adjEdgeIndex;
		}
	} else if (getCurrentFace(twin13) == faceIndex) {           // Case 3?
		// The extension fills a notch between vertices v2, v3, v1.
		index_t successor = getNextEdge(adjEdgeIndex);
		d_edgeTable[twin13].d_nextEdge = successor;
		d_edgeTable[successor].d_prevEdge = twin13;
		d_edgeTable[twin13].d_twinEdge = twin21;
		if (isDefined(twin21)) {
			d_edgeTable[twin21].d_twinEdge = twin13;
		}

        if (! isolatedVertexCheckAndDelete(adjEdgeIndex, successor)) {
          std::cerr << "Polyhedron<T>::extendFace(" << v1Index << ", "
                    << adjEdgeIndex << "): " 
                    << "Undefined edge encounted while traversing face "
                    << faceIndex << std::endl;
          return false;
        }
	} else {                                                    // Case 1.
		// The extension is open, only sharing one edge (adjEdgeIndex)
		// with the coplanar face.
		index_t nextEdgeIndex = edges();
		index_t successor = getNextEdge(adjEdgeIndex);
		d_edgeTable[adjEdgeIndex].d_nextEdge = nextEdgeIndex;
		d_edgeTable[adjEdgeIndex].d_twinEdge = twin13;
		d_edgeTable[successor].d_prevEdge = nextEdgeIndex;
		if (isDefined(twin13)) {
			d_edgeTable[twin13].d_twinEdge = adjEdgeIndex;
		}
		if (isDefined(twin21)) {
			d_edgeTable[twin21].d_twinEdge = nextEdgeIndex;
		}
		
		HalfEdge edge12(v1Index, adjEdgeIndex, successor, twin21, faceIndex);
		d_edgeTable.push_back(edge12);
		if (! isDefined(d_vertexTable[v1Index].d_edge)) 
			d_vertexTable[v1Index].d_edge = nextEdgeIndex;
	}
	
	return true;
}


template<typename T>
void Polyhedron<T>::computeNormalVectors() {
	for(index_t faceIndex = 0; faceIndex < faces(); ++faceIndex) {
		computeNormalVector(faceIndex);
	}
}


// Computes the normal vector using Newell's method
template<typename T>
void Polyhedron<T>::computeNormalVector(index_t faceIndex) {
	Face* face = &d_faceTable[faceIndex];
	FaceNormal &normal = face->d_normalVector;
	
	index_t edge = face->d_outer;
	index_t next = d_edgeTable[edge].d_nextEdge;
	normal.clear();
	
	do {
		FaceVertex& p = d_vertexTable[getCurrentVertex(edge)].d_v;
		FaceVertex& n = d_vertexTable[getNextVertex(edge)].d_v;
					
		normal[0] += (p[2] + n[2])*(p[1] - n[1]);
		normal[1] += (p[0] + n[0])*(p[2] - n[2]);
		normal[2] += (p[1] + n[1])*(p[0] - n[0]);
		
		edge = next;
		next = d_edgeTable[edge].d_nextEdge;
	} while (edge != face->d_outer);
	
	// normalize the normal vector only if S is double or float. Otherwise, divide
	// each integral component by the greatest common divisor. This ensures that the
	// contents of normal represent a three-dimensional vector that is perpendicular
	// to the face using the same numeric type as the vertices.
	normal.reduce();
}


// deleteFace marks the indicated face for removal, and functionally deletes
// the face from the table, including its edges, and any vertices that exclusively
// belong to deleted faces. 

template<typename T>
bool Polyhedron<T>::deleteFace2(index_t faceIndex) {
	assert(0 <= faceIndex && faceIndex < faces());
	
	// Identify the edges that delineate the selected face, and store them in the
	// list edgeLoop.
	index_t initialEdge = d_faceTable[faceIndex].d_outer;
	index_t edge = initialEdge;

    // The sets deadEdges and deadVertices are used to accumulate the indices
    // of edges and vertices to be deleted in decreasing order. (The default,
    // which uses std::less<index_t>, constructs the sets in increasing order.)
    IntSet deadEdges;
    IntSet deadVertices;

	do {
		if (! isDefined(edge)) {
			std::cerr << "Warning[deleteFace2(" << faceIndex << ")]: "
					  << "Attempted to delete an incomplete face."
					  << std::endl;
			return false;
		}
		
		deadEdges.insert(edge);

        // Ensure that each vertex in the loop, references an incident edge outside of
        // the loop. If this is not possible, then the vertex can be deleted.

        index_t vertexIndex = getCurrentVertex(edge);
		if (isDefined(vertexIndex)) {
			if (d_vertexTable[vertexIndex].d_edge == edge) {
				index_t altEdge;
				bool status = getCollocalEdge(edge, &altEdge);
				if (status) {
					d_vertexTable[vertexIndex].d_edge = altEdge;
				} else {
					deadVertices.insert(vertexIndex);
				}
			}
		}

		edge = getNextEdge(edge);	
	} while(edge != initialEdge);
	
	assert(edge == initialEdge);
	do {
		index_t vertexIndex = getCurrentVertex(edge);
		if (isDefined(vertexIndex)) {
			if (d_vertexTable[vertexIndex].d_edge == edge) {
				index_t altEdge;
				bool status = getCollocalEdge(edge, &altEdge);
				if (status) {
					d_vertexTable[vertexIndex].d_edge = altEdge;
				} else {
					deleteVertex2(vertexIndex);
				}
			}
		}
		edge = getNextEdge(edge);
	} while (edge != initialEdge);
	
	
    IntSet::iterator pos;
	for(pos = deadEdges.begin(); pos != deadEdges.end(); ++pos) {
		bool status = Polyhedron<T>::deleteEdge2(*pos);
		if (! status) {
			std::cerr << "Polyhedron<T>::deleteFace2(" << faceIndex << ")"
					<< " error encounterd during call to deleteEdge2(" << *pos << ")."
					<< std::endl;
			return false;
		}
	}

    for(pos = deadVertices.begin(); pos != deadVertices.end(); ++pos) {
		bool status = Polyhedron<T>::deleteVertex2(*pos);
		if (! status) {
			std::cerr << "Polyhedron<T>::deleteFace2(" << faceIndex << ")"
					<< " error encounterd during call to deleteVertex2(" << *pos << ")."
					<< std::endl;
			return false;
		}
	}
	
	index_t lastFace = faces() - 1;
	d_faceTable[faceIndex] = d_faceTable[lastFace];
	d_faceTable.pop_back();
	
	initialEdge = d_faceTable[faceIndex].d_outer;
	edge = initialEdge;
	do {
		if (! isDefined(edge)) {
			std::cerr << "Warning[deleteFace2(" << faceIndex << ")]: "
					  << "relocated an incomplete face."
					  << std::endl;
			return false;
		}
		d_edgeTable[edge].d_face = faceIndex;
		edge = getNextEdge(edge);
	} while (edge != initialEdge);
	return true;	
}

// getCollocalEdge searches for an edge different from edgeIndex that shares its vertex.
// If found, the index of the discovered collocal edge is returned via the pointer in
// the second argument, and the function returns the value true. If no edge can be
// found, then false is returned. Note that three search strategies are employed in
// succession.
template<typename T>
bool Polyhedron<T>::getCollocalEdge(index_t edgeIndex, index_t *collocalEdge) const {
	index_t vertexIndex = d_edgeTable[edgeIndex].d_vertex;
	*collocalEdge = getUndefined();
	
	if (isDefined(vertexIndex)) {
		*collocalEdge = getClockwiseEdge(edgeIndex); // returns an edge that is collocal with edgeIndex
												     // or UNDEFINED.
		if (! isDefined(*collocalEdge) || (! d_edgeTable[*collocalEdge].d_isPresent)) {
			*collocalEdge = getCounterClockwiseEdge(edgeIndex);
		}
		
		if (! isDefined(*collocalEdge) || (! d_edgeTable[*collocalEdge].d_isPresent)) {
			for(index_t edge = 0; edge < edges(); ++edge) {
				if (d_edgeTable[edge].d_vertex == vertexIndex && 
					edge != edgeIndex &&
					d_edgeTable[edge].d_isPresent) {
						*collocalEdge = edge;
						break;
					}
			}	
		}
	}
	return isDefined(*collocalEdge) && d_edgeTable[*collocalEdge].d_isPresent;
}

// delete the indicated edge and remove it from the d_edgeTable by overwriting
// d_edgeTable[edgeIndex] with the last item in the table, and then popping the
// d_edgeTable. Before calling this routine, all references to the edge to be
// removed should have been removed from the d_faceTable and d_vertexTable.
// Since the index of the last edge in the table has been changed from edges()
// to edgeIndex, all references to the last edge in the d_faceTable, d_edgeTable,
// and d_vertexTable should be updated.

template<typename T>
bool Polyhedron<T>::deleteEdge2(index_t edgeIndex) {
	// Delete the indicated edge by overwriting
	index_t lastEdgeIndex = edges() - 1;
	d_edgeTable[edgeIndex] = d_edgeTable[lastEdgeIndex];
	d_edgeTable.pop_back();
	
	// Update the face table, if needed.
	index_t faceIndex = d_edgeTable[edgeIndex].d_face;
	if (d_faceTable[faceIndex].d_outer == lastEdgeIndex) {
		d_faceTable[faceIndex].d_outer = edgeIndex;
	}
	
	// Update the vertex table, if needed.
	index_t vertexIndex = d_edgeTable[edgeIndex].d_vertex;
	if (d_vertexTable[vertexIndex].d_edge == lastEdgeIndex) {
		d_vertexTable[vertexIndex].d_edge = edgeIndex;
	}
	
	// Update the edge table:

	index_t nextEdge = getNextEdge(edgeIndex);
	if (getPrevEdge(nextEdge) == lastEdgeIndex) {
		d_edgeTable[nextEdge].d_prevEdge = edgeIndex;
	} else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
			<< "d_edgeTable[" << edgeIndex << "].d_nextEdge = "
			<< d_edgeTable[edgeIndex].d_nextEdge << " and "
			<< "d_edgeTable[" << nextEdge << "].d_prevEdge = "
			<< d_edgeTable[nextEdge].d_prevEdge;
		return false;
	}
	
	index_t twinEdge = getTwinEdge(edgeIndex);
	if (getTwinEdge(twinEdge) == lastEdgeIndex) {
		d_edgeTable[twinEdge].d_twinEdge = edgeIndex;
	} else if (isDefined(getTwinEdge(twinEdge))) {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
			<< "d_edgeTable[" << edgeIndex << "].d_twinEdge = "
			<< d_edgeTable[edgeIndex].d_twinEdge << " and "
			<< "d_edgeTable[" << twinEdge << "].d_twinEdge = "
			<< d_edgeTable[twinEdge].d_twinEdge;
		return false;		
	}
	
	index_t prevEdge = getPrevEdge(edgeIndex);
	if (getNextEdge(prevEdge) == lastEdgeIndex) {
		d_edgeTable[prevEdge].d_nextEdge = edgeIndex;
	} else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
			<< "d_edgeTable[" << edgeIndex << "].d_prevEdge = "
			<< d_edgeTable[edgeIndex].d_prevEdge << " and "
			<< "d_edgeTable[" << prevEdge << "].d_nextEdge = "
			<< d_edgeTable[prevEdge].d_nextEdge;
		return false;
	}
	
	return true;
}

// deleteVertex2 deletes the indicated vertex from d_vertexTable by replacing 
// entry vIndex with the last entry in the table. Before this procedure is
// called, all edges that reference vertex vIndex should have been deleted.
// All edges that radiate from the relocated vertex (those that are collocal),
// need to be updated with the new index.
template<typename T>
bool Polyhedron<T>::deleteVertex2(index_t vIndex) {
	index_t lastVertex = vertices() - 1;
	d_vertexTable[vIndex] = d_vertexTable[lastVertex];
	d_vertexTable.pop_back();
	
	return Polyhedron<T>::updateEdgesOfVertex(vIndex, lastVertex);	
}

// Function updateEdgesOfVertex replaces each reference to vertex vOld to the new value
// vNew in the entries in d_edgeTable. It is assumed that all edges that eminate from vertex
// vOld are topologically accessible by pointer algebra.

template<typename T>
bool Polyhedron<T>::updateEdgesOfVertex(index_t vNew, index_t vOld) {
  index_t initialEdge = d_vertexTable[vOld].d_edge;
  assert(0 <= vNew < vertices());
  assert(0 <= vOld < vertices());


  if (isDefined(initialEdge)) {		
    enum State {CLOCKWISE, COUNTERCLOCKWISE, SEQUENTIAL, FINISHED};
    State state = CLOCKWISE;
    index_t edge = initialEdge;
    
    while (state != FINISHED) {
      switch (state) {
      case CLOCKWISE:
        if (isPresentEdge(edge)) {
          d_edgeTable[edge].d_vertex = vNew;
          edge = getClockwiseEdge(edge);

          if (edge == initialEdge) {
            state = FINISHED;
          }
        } else {
          state = COUNTERCLOCKWISE;
          edge = getCounterClockwiseEdge(initialEdge);
        }
        break;

      case COUNTERCLOCKWISE:
        if (isPresentEdge(edge)) {
          d_edgeTable[edge].d_vertex = vNew;
          edge = getCounterClockwiseEdge(edge);

          if (edge == initialEdge) {
            state = FINISHED;                                   // This seems unlikely.
          }
        } else {
          state = SEQUENTIAL;
        }
        break;

      case SEQUENTIAL:
        // Replacement was not found by edge traversal. This likely means
        // that the last face adjacent to this vertex is being marked
        // for deletion. We can confirm this by traversing the
        // edge table.
        for(edge = 0; edge < edges(); edge++) {
          if (d_edgeTable[edge].d_vertex == vOld && isPresentEdge(edge)) {
            d_edgeTable[edge].d_vertex = vNew;
          }
        }
        state = FINISHED;
        break;
      }
    }
    return true;
  }
  return false;
}

// After all faces have been deleted purgeFaces may be called to reindex the remaining
// polyhedral elements, and remove the marked faces, edges, and vertices from the tables 
template<typename T>
bool Polyhedron<T>::deleteFace(index_t faceIndex) {
	assert(0 <= faceIndex && faceIndex < faces());
	
	if (d_faceTable[faceIndex].d_isPresent == false) {
		std::cerr << 
			"Polyhedron<T>::deleteFace(faceIndex = " << faceIndex << 
			"): WARNING: attempt to delete a deleted face."	<< std::endl;
		return false;
	}
	
	d_faceTable[faceIndex].d_isPresent = false;  // marks the face for deletion
	d_deadFaces.insert(faceIndex);

	index_t initialEdge = d_faceTable[faceIndex].d_outer;
	if (isDefined(initialEdge)) {
		index_t edge = initialEdge;
		do {
          // Delete each edge of the face marked for deletion.
			bool status = Polyhedron<T>::deleteEdge(edge);
            if (! status) {
              return false;
            }

			edge = getNextEdge(edge); // walk-around the face, counter-clockwise.
			if (! isDefined(edge)) {
				std::cerr << "Warning[deleteFace(" << faceIndex << ")]: "
						  << "Attempted to delete an incomplete face."
						  << std::endl;
				return false;
			}		
		} while (edge != initialEdge);
	}
	return true;
}

// deleteEdge marks the indicated edge for deletion
template<typename T>
bool Polyhedron<T>::deleteEdge(index_t edgeIndex) {
	assert(0 <= edgeIndex && edgeIndex < edges());
	
	if (d_edgeTable[edgeIndex].d_isPresent == false) {
		std::cerr << 
			"Polyhedron<T>::deleteEdge(edgeIndex = " << edgeIndex << 
			"): WARNING: attempt to delete a deleted edge."	<< std::endl;
		return false;
	}
	
	d_edgeTable[edgeIndex].d_isPresent = false;
	d_deadEdges.insert(edgeIndex);

	// remove any reference to edge by its twin.
	index_t twinEdge = getTwinEdge(edgeIndex);
	if (isDefined(twinEdge)) {
		if (d_edgeTable[twinEdge].d_twinEdge == edgeIndex) {
			d_edgeTable[twinEdge].d_twinEdge = UNDEFINED; 
		} else {
			std::cerr <<
				"Polyhedron<T>::deleteEdge(edgeIndex = " << edgeIndex << 
				"): ERROR: twin inconsistency with edge "
				<< twinEdge << std::endl;
            return false;
		}
	}		
	
    // Modify d_faceTable references to edgeIndex:
	index_t faceIndex = d_edgeTable[edgeIndex].d_face;
	if (d_faceTable[faceIndex].d_outer == edgeIndex && 
		d_faceTable[faceIndex].d_isPresent) {
		// If possible, replace d_faceTable[faceIndex].d_outer with
		// an edge that is retained (or present).
		
		index_t edge = getNextEdge(edgeIndex);
		do {
			if (d_edgeTable[edge].d_face == faceIndex &&
				d_edgeTable[edge].d_isPresent) {
					d_faceTable[faceIndex].d_outer = edge;
					break;
			}
			edge = getNextEdge(edge);
		} while (edge != edgeIndex && isDefined(edge));

        if (d_faceTable[faceIndex].d_outer == edgeIndex) {
          // Could not update d_faceTable[faceIndex].d_outer with an active edge.
          d_faceTable[faceIndex].d_outer = UNDEFINED;
        }
	}
	
    // Modify d_vertexTable references to edgeIndex:
	index_t vertexIndex = d_edgeTable[edgeIndex].d_vertex;
	if (d_vertexTable[vertexIndex].d_edge == edgeIndex) {
		// If possible, replace d_vertexTable[d_edgeTable[edge].d_vertex].d_edge
		// with an active edge that originates from vertexIndex.
		index_t edge = getClockwiseEdge(edgeIndex);
		enum State {CLOCKWISE, COUNTERCLOCKWISE, SEQUENTIAL, SUCCESS, FAILURE, FINISHED};
		State state = CLOCKWISE;
		while (state != FINISHED) {
			switch (state) {
				case CLOCKWISE:
					if (! isDefined(edge)) {
						state = COUNTERCLOCKWISE;
						edge = getCounterClockwiseEdge(edgeIndex);
					} else if (edge == edgeIndex) {
						state = SEQUENTIAL;
					} else if (d_edgeTable[edge].d_isPresent) {
						state = SUCCESS;	
					} else {
						edge = getClockwiseEdge(edge);
					}
					break;
				case COUNTERCLOCKWISE:
					if (! isDefined(edge)) {
						state = SEQUENTIAL;
						edge = edgeIndex;
					} else if (edge == edgeIndex) { // this should not occur!
						state = SEQUENTIAL;
					} else if (d_edgeTable[edge].d_isPresent) {
						state = SUCCESS;
					} else {
						edge = getCounterClockwiseEdge(edge);
					}
					break;
				case SEQUENTIAL:
					// Replacement was not found by edge traversal. This likely means
					// that the last face adjacent to this vertex is being marked
					// for deletion. We can confirm this by traversing the
					// edge table.
					for(edge = 0; edge < edges(); edge++) {
						if (d_edgeTable[edge].d_vertex == vertexIndex && d_edgeTable[edge].d_isPresent) {
							state = SUCCESS;
							break;
						}
					}
					if (state != SUCCESS) {
						state = FAILURE;
					}
					break;
					
				case SUCCESS:
					d_vertexTable[vertexIndex].d_edge = edge;
					state = FINISHED;
					break;
				
				case FAILURE:
					Polyhedron<T>::deleteVertex(vertexIndex);
					state = FINISHED;
					break;
			}
		}
	}
	return true;
}

template<typename T>
bool Polyhedron<T>::deleteVertex(index_t vertexIndex) {
	assert(0 <= vertexIndex && vertexIndex < vertices());
	d_vertexTable[vertexIndex].d_edge = UNDEFINED;
	d_vertexTable[vertexIndex].d_isPresent = false;
    d_deadVertices.insert(vertexIndex);
	return true;
}
	
// purgeFaces removes all faces that are marked for deletion, i.e., those
// faces that contain a value of d_isPresent == false. In addition, all halfEdges
// and vertices that are dependent on a face that is removed are deleted.
template<typename T>
bool Polyhedron<T>::purgeFaces(void) {
  IntSet::iterator pos;
  for (pos = d_deadFaces.begin(); pos != d_deadFaces.end(); ++pos) {
    index_t lastFaceIndex = d_faceTable.size() - 1;
    if (*pos != lastFaceIndex) {
      // fill position *pos in the table with the last face.
      d_faceTable[*pos] = d_faceTable[lastFaceIndex];
      // update cross references from edgeTable
      index_t initialEdge = d_faceTable[*pos].d_outer;
      index_t edge = initialEdge;
      do {
		if (! isDefined(edge)) {
          std::cerr << "Polyhedron<T>::purgeFaces: "
                    << "d_faceTable[" << lastFaceIndex << "] is incomplete."
                    << std::endl;
          return false;
		}

		d_edgeTable[edge].d_face = *pos;
		edge = getNextEdge(edge);
      } while (edge != initialEdge);
    }
    d_faceTable.pop_back();
  }
  d_deadFaces.clear();
  return true;
}

// Function purgeVertices removes from d_vertexTable all records of vertices that have
// been marked for deletion. A vertex is marked for deletion by invoking function
// deleteVertex. The latter inserts the index value of the marked vertex into the
// buffer d_deadVertices.
template<typename T>
bool Polyhedron<T>::purgeVertices(void) {
 IntSet::iterator pos;
  for (pos = d_deadVertices.begin(); pos != d_deadVertices.end(); ++pos) {
    index_t lastIndex = d_vertexTable.size() - 1; // the last (active) face in the table.
    if (*pos != lastIndex) {
      // Move the last face (lastIndex) into position *pos in the d_vertexTable.

      // First copy the last face entries into position *pos. This purges vertex *pos.
      d_vertexTable[*pos] = d_vertexTable[lastIndex];

      // Then change all references to vertex lastIndex that appear in the d_edgeTable
      // to vertex *pos.
      bool status = updateEdgesOfVertex(*pos, lastIndex);
      if (! status) {
        std::cerr << "Polyhedron<T>:purgeVertices: "
                  << "unable to update edges of vertex " << lastIndex
                  << std::endl;
        return false;
      }
    }

    // Finally delete the last (active) entry, since it is now redundant.
    d_vertexTable.pop_back();
  }

  // Once all vertices have been purged, clear the d_deadVertices buffer.
  d_deadVertices.clear();

  // Recompute the bounding box of the polyhedron.
  computeVertexExtrema();

  return true;
}

template<typename T>
void Polyhedron<T>::computeVertexExtrema(void) {
  if (d_vertexTable.size() > 0) {
	typename std::vector<Vertex>::iterator v;
	d_vMin = d_vertexTable[0].d_v;
	d_vMax = d_vertexTable[0].d_v;
	for (index_t i = 1; i < d_vertexTable.size(); ++i) {
	  for(size_t j = 0; j < D; j++) {
		d_vMin[j] = min(d_vMin[j], d_vertexTable[i].d_v[j]);
		d_vMax[j] = max(d_vMax[j], d_vertexTable[i].d_v[j]);
	  }
	}
  }
}

template<typename T>
bool Polyhedron<T>::purgeEdges(void) {
 IntSet::iterator pos;
  for (pos = d_deadEdges.begin(); pos != d_deadEdges.end(); ++pos) {
    index_t lastIndex = d_edgeTable.size() - 1;
    if (*pos != lastIndex) {
      // fill position *pos in the table with the last face.
      d_edgeTable[*pos] = d_edgeTable[lastIndex];

      // Update the face table, if needed.
      index_t faceIndex = d_edgeTable[*pos].d_face;
      if (d_faceTable[faceIndex].d_outer == lastIndex) {
		d_faceTable[faceIndex].d_outer = *pos;
      }
	
      // Update the vertex table, if needed.
      index_t vertexIndex = d_edgeTable[*pos].d_vertex;
      if (d_vertexTable[vertexIndex].d_edge == lastIndex) {
		d_vertexTable[vertexIndex].d_edge = *pos;
      }
	
      // Update the edge table:

      index_t nextEdge = getNextEdge(*pos);
      if (getPrevEdge(nextEdge) == lastIndex) {
		d_edgeTable[nextEdge].d_prevEdge = *pos;
      } else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
                  << "d_edgeTable[" << *pos << "].d_nextEdge = "
                  << d_edgeTable[*pos].d_nextEdge << " and "
                  << "d_edgeTable[" << nextEdge << "].d_prevEdge = "
                  << d_edgeTable[nextEdge].d_prevEdge;
		return false;
      }
	
      index_t twinEdge     = getTwinEdge(*pos);
      index_t twinTwinEdge = getTwinEdge(twinEdge);
      if (isDefined(twinTwinEdge)) {
        if (twinTwinEdge == lastIndex) {
          d_edgeTable[twinEdge].d_twinEdge = *pos;
        } else {
          std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
                    << "d_edgeTable[" << *pos << "].d_twinEdge = "
                    << d_edgeTable[*pos].d_twinEdge << " and "
                    << "d_edgeTable[" << twinEdge << "].d_twinEdge = "
                    << d_edgeTable[twinEdge].d_twinEdge;
          return false;		
        }
      }
	
      index_t prevEdge = getPrevEdge(*pos);
      if (getNextEdge(prevEdge) == lastIndex) {
		d_edgeTable[prevEdge].d_nextEdge = *pos;
      } else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
                  << "d_edgeTable[" << *pos << "].d_prevEdge = "
                  << d_edgeTable[*pos].d_prevEdge << " and "
                  << "d_edgeTable[" << prevEdge << "].d_nextEdge = "
                  << d_edgeTable[prevEdge].d_nextEdge;
		return false;
      }
    }
    d_edgeTable.pop_back();
  }

  d_deadEdges.clear();
  return true;
}


template<typename T>
bool Polyhedron<T>::purgeTables(void) {
  bool status = purgeFaces();
  if (! status) {
    std::cerr << "Polyhedron<T>::purgeTables: ERROR: "
              << "Failed to purgeFaces."
              << std::endl;
    return false;
  }
 
  status = purgeVertices();
  if (! status) {
    std::cerr << "Polyhedron<T>::purgeTables: ERROR: "
              << "Failed to purgeVertices."
              << std::endl;
    return false;
  }
 
  status = purgeEdges();
  if (! status) {
    std::cerr << "Polyhedron<T>::purgeTables: ERROR: "
              << "Failed to purgeEdges."
              << std::endl;
    return false;
  }

  return true;
}

template<typename T>
void Polyhedron<T>::describe(std::ostream &os) const {
	const int formatWidth = 12;
	index_t nVertices = vertices();
	index_t nEdges    = edges();
	index_t nFaces    = faces();
	
	os << "Vertex Table (" << nVertices << ")" << std::endl;
	os 	<< std::setw(formatWidth) << "V " 
		<< std::setw(formatWidth*2) << "Position " 
		<< std::setw(formatWidth - 3) << "Edge "
		<< std::setw(formatWidth) << "Present (T/F) " << std::endl;
		
	for(index_t i = 0; i < nVertices; ++i) {
		const Vertex &v = d_vertexTable[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << std::setprecision(8) << v.d_v 
			<< std::setw(formatWidth) << printIndex(v.d_edge)
			<< std::setw(formatWidth) << v.d_isPresent << std::endl;
	}
	
	os << std::endl<< "Edge Table (" << nEdges << ")" << std::endl;
	os 	<< std::setw(formatWidth)     << "       E " 
		<< std::setw(formatWidth+3)   << "  Vertex " 
		<< std::setw(formatWidth+1)   << "NextEdge " 
		<< std::setw(formatWidth)     << "PrevEdge " 
		<< std::setw(formatWidth)     << "TwinEdge " 
		<< std::setw(formatWidth-3)   << "    Face " 
		<< std::setw(formatWidth+4)   << "Present (T/F) " << std::endl;
	for(index_t i = 0; i < nEdges; ++i) {
		const HalfEdge &e = d_edgeTable[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << printIndex(e.d_vertex)
			<< std::setw(formatWidth) << printIndex(e.d_nextEdge)
		   	<< std::setw(formatWidth) << printIndex(e.d_prevEdge)
			<< std::setw(formatWidth) << printIndex(e.d_twinEdge)
			<< std::setw(formatWidth) << printIndex(e.d_face)
			<< std::setw(formatWidth) << e.d_isPresent << std::endl;
	}
	
	os 	<< std::endl << "Face Table (" << nFaces << ") " << std::endl;
	os 	<< std::setw(formatWidth) << "F " 
		<< std::setw(formatWidth+4) << "OuterEdge " 
		<< std::setw(formatWidth+4)   << "NormalVector "
		<< std::setw(formatWidth+4)   << "Present (T/F) " << std::endl;
	
	for(index_t i = 0; i < nFaces; ++i) {
		const Face &f = d_faceTable[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << printIndex(f.d_outer)
			<< std::setw(formatWidth) << f.d_normalVector 
			<< std::setw(formatWidth) << f.d_isPresent<< std::endl;
	}
	
	os 	<< "Volume = " << this->volume() 
		<< "; Surface Area = " << this->surfaceArea() << std::endl;
	
}

// isConvex() returns true if the polyhedron is convex; false, if it is nonconvex.

template<typename T>
bool Polyhedron<T>::isConvex() const {

	typename std::vector<typename Polyhedron<T>::Face>::const_iterator facePos;
	for(facePos = d_faceTable.begin(); facePos != d_faceTable.end(); ++facePos) {
		if (facePos->d_isPresent) {
			FaceVertex faceCenter; // The center of the current face.
			int vertexCount;		// Number of vertices that define the current face.

			// Compute the face center.
			faceCenter.clear();  // Initialize accumulators.
			vertexCount = 0;

			index_t initialEdge = facePos->d_outer;  // Identify the initial edge.
			index_t edge = initialEdge;
			do {	// for each edge that delimits the current face.
				faceCenter += d_vertexTable[getCurrentVertex(edge)].d_v;
				vertexCount++;
				edge = getNextEdge(edge);
				if (edge == UNDEFINED) {
					std::cerr << "WARNING[Polyhedron<T>::isConvex()]: "
						<< "Polyhedron has an incomplete face." << std::endl;
					return false;
				}	
			} while (edge != initialEdge);
			faceCenter /= vertexCount;

			typename std::vector<typename Polyhedron<T>::Vertex>::const_iterator vertexPos;
			for(vertexPos = d_vertexTable.begin(); vertexPos != d_vertexTable.end(); ++vertexPos) {
				FaceVertex displacement = vertexPos->d_v - faceCenter;
				T ip = facePos->d_normalVector * displacement;  // inner product
				if (ip  > 0) {
					std::cout << "Nonconvex relation: " 
						<< facePos->d_normalVector << "* ("
						<< vertexPos->d_v << " - " 
						<< faceCenter << ").normalize()" << " =  "
						<< ip << std::endl;
					return false;
				}
			}
		}
	}
	return true;
}

template<typename T>
double Polyhedron<T>::faceArea(index_t faceIndex) const {
	FaceVertex edgeProd;
	index_t initialEdge = d_faceTable[faceIndex].d_outer;
	index_t edge = initialEdge;
	do {
		edgeProd += crossProduct(d_vertexTable[getCurrentVertex(edge)].d_v, 
								 d_vertexTable[getNextVertex(edge)].d_v);
		edge = getNextEdge(edge);
	} while (edge != initialEdge);
	
	return (edgeProd*d_faceTable[faceIndex].d_normalVector)/2.0;	// Inner product		
}

template<typename T>
double Polyhedron<T>::surfaceArea() const {
	T area = 0;
	for (index_t i = 0; i < faces(); i++) {
		if (d_faceTable[i].d_isPresent) {
			area += faceArea(i);
		}
	}
	return area;
}

template<typename T>
index_t Polyhedron<T>::countInactiveFaces() const {
  index_t count = 0;
  for(index_t i = 0; i < d_faceTable.size(); ++i) {
    if (! d_faceTable[i].d_isPresent) count++;
  }
  return count;
}

template<typename T>
index_t Polyhedron<T>::countInactiveEdges() const {
  index_t count = 0;
  for(index_t i = 0; i < d_edgeTable.size(); ++i) {
    if (! d_edgeTable[i].d_isPresent) count++;
  }
  return count;
}

template<typename T>
index_t Polyhedron<T>::countInactiveVertices() const {
  index_t count = 0;
  for(index_t i = 0; i < d_vertexTable.size(); ++i) {
    if (! d_vertexTable[i].d_isPresent) count++;
  }
  return count;
}

template<typename T>
bool Polyhedron<T>::checkTables(bool verbose) const {
  bool status = true;

  index_t xFaces = countInactiveFaces();
  if (xFaces > 0) {
    std::cerr << "Warning: " << xFaces << " are inactive." << std::endl;
    status = false;
  }

  index_t xEdges = countInactiveEdges();
  if (xEdges > 0) {
    std::cerr << "Warning: " << xEdges << " are inactive." << std::endl;
    status = false;
  }

  index_t xVertices = countInactiveVertices();
  if (xVertices > 0) {
    std::cerr << "Warning: " << xVertices << " are inactive." << std::endl;
    status = false;
  }

  index_t dFaces = d_deadFaces.size();
  if (dFaces != xFaces) {
    std::cerr << "Warning: number of deadFaces (" << dFaces << ") "
              << "disagrees with number of inactive faces (" << xFaces << ")."
              << std::endl;
    status = false;
  }

  index_t dEdges = d_deadEdges.size();
  if (dEdges != xEdges) {
    std::cerr << "Warning: number of deadEdges (" << dEdges << ") "
              << "disagrees with number of inactive edges (" << xEdges << ")."
              << std::endl;
    status = false;
  }

  index_t dVertices = d_deadVertices.size();
  if (dVertices != xVertices) {
    std::cerr << "Warning: number of deadVertices (" << dVertices << ") "
              << "disagrees with number of inactive vertices (" << xVertices << ")."
              << std::endl;
    status = false;
  }

  return status;
}
    
// Function eulerTest tests the topological consistency of the polyhedron defined
// by the three tables: d_faceTable, d_edgeTable, and d_vertexTable. If the argument
// verbose is true, then a detailed error report is printed on std::cerr. If no
// topological errors are detected, then eulerTest returns true, otherwise false.

template<typename T>
bool Polyhedron<T>::eulerTest(bool verbose) const {
  std::valarray<index_t> eFreq(static_cast<index_t>(0), d_edgeTable.size());
  std::valarray<index_t> vFreq(static_cast<index_t>(0), d_vertexTable.size());
  index_t faceCount = 0;
  index_t edgeCount = 0;
  index_t vertexCount = 0;
  index_t errorCount = 0;
	
  // for each active face
  for(index_t i = 0; i < d_faceTable.size(); i++) {
    if (d_faceTable[i].d_isPresent) {
      faceCount++;
			
      index_t firstEdge = d_faceTable[i].d_outer;
      index_t edge      = firstEdge;
      index_t twin;
      do {

        if (! isDefined(edge) ) {
          std::cerr << "eulerTest: ERROR : face " << i << " is incomplete!"
                    << " Cannot continue test!"
                    << std::endl;
          return false;
        }

        if (! d_edgeTable[edge].d_isPresent) {
          if (verbose) {
            std::cerr << "eulerTest: Edge " << edge << ", required by face " 
                      << i << " is marked for deletion." 
                      << std::endl;
          }
          errorCount++;
        }

        if ( d_edgeTable[edge].d_face != i) {
          if (verbose) {
            std::cerr <<  "eulerTest: Edge " << edge << ", belongs to face " 
                      << i << " but references face " << d_edgeTable[edge].d_face
                      << std::endl;
          }
          errorCount++;
        }
				
        // test twin consistency
        twin = getTwinEdge(edge);
        if (! isDefined(twin)) {
          if (verbose) {
            std::cerr << "eulerTest: Edge " << edge << " has no twin!" 
                      << std::endl;
          }
          errorCount++;
        } else if (getTwinEdge(twin) != edge) {
          if (verbose) 
            std::cerr << "eulerTest:WARNING: The twin of edge " << edge 
                      << " has edge " << d_edgeTable[twin].d_twinEdge 
                      << " as its twin!" << std::endl;
          errorCount++;
        }

        // test next & prev consistency
        index_t nextEdge = getNextEdge(edge);
        if (getPrevEdge(nextEdge) != edge) {
          if (verbose) {
            std::cerr << "eulerTest: Pointer inconsistency between edge " << edge
                      << " and its successor " << nextEdge 
                      << std::endl;
          }
          errorCount++;
        }

				
        eFreq[edge]++;                                  // indicate that edge was used.
        vFreq[d_edgeTable[edge].d_vertex]++;              // indicate that vertex was used
				
    		
        edge = nextEdge;
      } while (edge != firstEdge);
    }
  }
	
  // Report on edge usage.
  for(index_t i = 0; i < eFreq.size(); i++) {
    if (eFreq[i] == 1) {
      edgeCount++;
    } else if (eFreq[i] > 1) {
      if (verbose) {
        std::cerr << "eulerTest: WARNING: edge " << i << "was used " 
                  << eFreq[i] << " times."
                  << std::endl;
      }
      errorCount++;
    } else if (d_edgeTable[i].d_isPresent) {
      if (verbose) 
        std::cerr << "eulerTest:WARNING: Edge " 
                  << i << " is present but not used." << std::endl;
      errorCount++;
    }
  }
	
  if (edgeCount % 2 == 1) {
    if (verbose) 
      std::cerr << "eulerTest:WARNING: the number of half edges equals an odd integer." 
                << std::endl;
    errorCount++;
  }
	
  for(index_t i = 0; i < vFreq.size(); i++) {
    if (vFreq[i] > 0) {
      vertexCount++;
    }
  }

  int eulerInvariant = vertexCount - edgeCount/2 + faceCount;
  if (eulerInvariant != 2) errorCount++;
	
  if (verbose) 
    std::cout << "eulerTest: vertices = " << vertexCount 
              << ", half edges = " << edgeCount 
              << ", faces = " << faceCount << std::endl
              << "V - HE/2 + F = " << eulerInvariant << std::endl;
			
  return (errorCount == 0);		
}


template<typename T>
  double Polyhedron<T>::volume() const {
  double volume = 0;
  FaceVertex centroid;
  index_t vertexCount = 0;
	
  typename std::vector<Vertex>::const_iterator vertexPos;
  for(vertexPos = d_vertexTable.begin(); vertexPos != d_vertexTable.end(); ++vertexPos) {
    centroid += vertexPos->d_v;
    vertexCount++;
  }
  centroid /= vertexCount;
	
  typename std::vector<Face>::const_iterator facePos;
  for(facePos = d_faceTable.begin(); facePos != d_faceTable.end(); ++facePos) {
    if (facePos->d_isPresent) {
      FaceVertex edgeProd;
      FaceVertex faceCenter;
      index_t count = 0;
		
      edgeProd.clear();
      faceCenter.clear();
		
      index_t initialEdge = facePos->d_outer;
      index_t edge = initialEdge;
      do {
        count++;
        faceCenter += d_vertexTable[getCurrentVertex(edge)].d_v;
        edgeProd   += crossProduct(	d_vertexTable[getCurrentVertex(edge)].d_v, 
                                    d_vertexTable[getNextVertex(edge)].d_v);
        edge = getNextEdge(edge);
      } while (edge != initialEdge);
			
      faceCenter /= count;
			
      double sectionVolume = edgeProd*(faceCenter - centroid); // inner product	
      volume += sectionVolume;
    }
  }
  return volume/6.0; // Geometric factor (1/2)*(1/3)
}

template<typename T>
void Polyhedron<T>::display() const {
  typename std::vector<Face>::const_iterator facePos;
  for(facePos = d_faceTable.begin(); facePos != d_faceTable.end(); ++facePos) {
	if (facePos->d_isPresent) {
	  index_t edge = facePos->d_outer;

	  GLdouble nv[D];
	  for(int i = 0; i < D; i++) {
		nv[i] = static_cast<GLdouble>(facePos->d_normalVector[i]);
	  }
	  glBegin(GL_POLYGON);

	  glNormal3dv(nv);

	  do {
		GLdouble v[D];
		index_t vertex = d_edgeTable[edge].d_vertex;

		for(int i = 0; i < D; i++) {
		  v[i] = static_cast<GLdouble>(d_vertexTable[vertex].d_v[i]);
		}
		glVertex3dv(v);

		edge = getNextEdge(edge);
	  } while (edge != facePos->d_outer && isDefined(edge));

	  glEnd();
	}
  }
}


#ifdef COMMENT
template <typename T>
std::ostream<charT,traits>& operator<< (std::ostream &os, 
                                        Polyhedron<T>::HalfEdge &e) {
  std::ostringstream s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision);

  s << "HalfEdge(" 
    << e.d_face << ", " 
    << e.d_prevEdge << ", " 
    << e.d_nextEdge << ", " 
    << e.d_twinEdge << ", " 
     << e.d_vertex << ","
     << e.d_isPresent << ")";
  os << s.str();
  return os;
}
#endif
  
#endif // __POLYHEDRON_H__
