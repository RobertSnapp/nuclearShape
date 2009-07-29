#ifndef __POLYGON_H__
#define __POLYGON_H__
//	$Id$  

//
//  Polygon.h
//
//  Created by Robert Snapp on 2007-07-10.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include <algorithm>
#include <list>
#include <set>
#include <vector>
#include <iterator>
#include "ntuple.h"
#include <iostream>
#include <iomanip>
#include <valarray>
#include "env.h"

typedef double Point2d[2];
typedef float  Point2f[2];
typedef long   Point2i[2];
typedef size_t index_t;    // index_t will be used to refer to edge and vertex indices.

// The parameter epsilonFind is used to test for equality between floating point entities. Thus, two
// vertices are considered to be coincident if the lie within epsilonFind of one another.
double const epsilonFind = 1.0E-08;

// A polygon is an ordered set of two-dimensional vertices and edges.
// The typename T designates the type of each vertex component (e.g.,
// char, short, int, float, or double.)
template<typename T>
class Polygon {


public:
  // typedef double normalType;
  static int const dim = 2;
  typedef ntuple<double,dim>  EdgeNormal;
  typedef ntuple<T,dim>       PolygonVertex;
  typedef std::set<index_t, std::greater<index_t> >  IntSet;

  static const int SEARCH_FAILED = -1;
	
  Polygon() : initialEdge_(UNDEFINED) {} // default constructor
  Polygon(std::vector<PolygonVertex> const &points);

  virtual ~Polygon() {}  // So this class can serve as a base class.

  // copy constructor
  Polygon(Polygon<T> const &p) {
    edgeTable_   = p.edgeTable_;
    vertexTable_ = p.vertexTable_;
    initialEdge_ = p.initialEdge_;
  }
	
  // addEdge adds a new edge to the polygon model; the index of the new edge is returned through
  // the pointer in the third argument.
  bool addEdge(index_t vertex, index_t prevEdge, index_t nextEdge, index_t* edgeIndex);
  bool addEdge(PolygonVertex &vertex, index_t prevEdge, index_t nextEdge, index_t* edgeIndex);

  bool createTriangle(PolygonVertex &v0, PolygonVertex &v1, PolygonVertex &v2);
	
  // marks the indicated edge (along with its implicit vertices,
  // if they are unique) for deletion from the current model.
  // deleteEdge returns true if the face was deleted without detecting any
  // errors. If false is returned, then an error was detected, and the 
  // calling routine should either throw an exception or abort.
  bool deleteEdge(index_t edgeIndex);
	
  // marks the indicated vertex for deletion.
  bool deleteVertex(index_t vertexIndex);
	
  // removes all faces, halfedges, and vertices that are marked for deletion.
  // purgeFaces returns the number of faces that were purged.
  bool purgeEdges();
  bool purgeVertices();
  bool purgeTables();
		
	
  // Returns the number of edges in the current polytope.
  index_t edges() const    {return static_cast<index_t>(edgeTable_.size());  }
	
  // Returns the number of vertices in the current polytope.
  index_t vertices() const {return static_cast<index_t>(vertexTable_.size());}
	
  // isConvex() returns true of the polygon is convex.
  bool isConvex(bool verbose) const;
	
  // Returns the next available index for the vertex table.
  index_t nextVertexIndex() {return vertices();}
	
  // Returns the next available index for the edge table.
  index_t nextEdgeIndex() {return edges();}
	
  // Returns the last index in the edge table.
  index_t lastEdgeIndex() {return edges() - 1;}
	
  EdgeNormal &getEdgePerp(index_t index) {
    assert(0 <= index && index < edges());
    index_t nextEdge = getNextEdge(index);

    EdgeNormal delta = vertexTable_[edgeTable_[nextEdge].vertex_].v_ 
      - vertexTable_[edgeTable_[index].vertex_].v_;
    return perp(delta);
  }
  

  bool isVisible(index_t edge, const PolygonVertex &p) const{
    if (! isDefined(edge))
      return false;

    index_t next = getNextEdge(edge);
    if (! isDefined(next)) 
      return false;

    index_t v1Index = getCurrentVertex(edge);
    if (! isDefined(v1Index))
      return false;

     index_t v2Index = getCurrentVertex(next);
     if (! isDefined(v2Index))
      return false;

    PolygonVertex v2 = vertexTable_[v2Index].v_;
    PolygonVertex v1 = vertexTable_[v1Index].v_;
    PolygonVertex delta = v2 - v1;
    return (p - v1)*perp(delta) > 0;
  }
  	
  // Checks the two tables (edgeTable_, and vertexTable_) for the presence of
  // inactive elements. If no inactive elements are detected, then checkTables returns
  // true, otherwise false.
  bool checkTables() const;

  // Function countInactiveEdges returns the number of edges in the edgeTable_ that are not
  // active, i.e. have present_ set to false.
  int countInactiveEdges() const;

  // Function countInactiveVertices returns the number of vertices in the vertexTable_ that are not
  // active, i.e. have present_ set to false.
  int countInactiveVertices() const;

  // Function eulerTest tests the topological consistency of the polyhedron defined
  // by the three tables: faceTable_, edgeTable_, and vertexTable_. If the argument
  // verbose is true, then a detailed error report is printed on std::cerr. If no
  // topological errors are detected, then eulerTest returns true, otherwise false.
  bool topologyTest(bool verbose = false) const; 
	
  // Prints a description of the current polyhedron on the indicated stream.
  void describe(std::ostream &os) const;
	
  // Computes are returns the area of the face with index i.
  double area() const;
 
  // Computes and returns the surface are of the current polyhedron.
  double perimeter() const ;

  void drawPerimeter() const;
	
protected:
  class Edge;   // Nested class for representing each HalfEdge
  class Vertex; // Nested class for representing each Vertex

	
  // Serves as a null pointer.
  static const index_t UNDEFINED = static_cast<index_t>(-1);
  static const int poly_level = 2; // debug level
	
#if (DEBUG 	> 0)	
  static const bool debug = (DEBUG >= poly_level);
#else
  static const bool debug = false;
#endif
	
  index_t getUndefined() const {return UNDEFINED;}
	
  // Returns true if a halfedge already exists that extends from vertex p 
  // to vertex n, within precision epsilon. If true, then index returns the
  // location of this edge in the HalfEdge table.
  bool findEdge(PolygonVertex const &p, PolygonVertex const &n, index_t *index, 
                double epsilon = 0) const;
	
  // findEdge(int, int, int) returns true if there is an edge in the edge
  // table that extends from the initial vertex to the final vertex. In
  // this case, the index of the discovered edge is returned via the 
  // third argument. If no edge is found, false is returned.
  bool findEdge(index_t initialVertexIndex, index_t finalVertexIndex, index_t *index) const;
	
  // Returns true if a vertex v already exists in the vertex table, within 
  // precision epsilon. If true, then index returns the location of this 
  // vertex in the vertex table.
  bool findVertex(PolygonVertex const &v, index_t* index, double epsilon = 0) const;
		
  // Returns the index of the current Vertex for the indicated edge.
  index_t getCurrentVertex(index_t edge) const {
    return (isDefined(edge) ? edgeTable_[edge].vertex_ : getUndefined());
  }
	
  // Returns the edge index that follows the value edge.
  index_t getNextEdge(index_t edge) const {
    return (isDefined(edge) ? edgeTable_[edge].nextEdge_ : getUndefined());
  }
	
 	
  // Returns the index of the next Vertex for the indicated edge.
  index_t getNextVertex(index_t edge) const {
    return getCurrentVertex(getNextEdge(edge));
  }
	
  // Returns the previous edge
  index_t getPrevEdge(index_t edge) const {
    return (isDefined(edge) ? edgeTable_[edge].prevEdge_ : getUndefined());
  }
		
  // Returns the edge index that precedes the value edge.
  index_t getPrevVertex(index_t edge) const {
    return getCurrentVertex(getPrevEdge(edge));
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

  bool isPresentEdge(index_t index) const {
    return (0 <= index && index < edges() && edgeTable_[index].present_);
  }
	
  index_t edgeOfVertex(index_t vertexIndex) {
    return (isDefined(vertexIndex) ? vertexTable_[vertexIndex].edge_ : getUndefined());
  }
	
  bool vertexIsIsolated(index_t vertexIndex) {
    return edgeOfVertex(vertexIndex) == UNDEFINED; 
  }
	
		
  // Clear the polytope.
  void clear() {
    edgeTable_.clear();
    vertexTable_.clear();
    initialEdge_ = UNDEFINED;
  }
	
  // ================
  // = Data members =
  // ================
  index_t               initialEdge_;
  std::vector<Edge>     edgeTable_;
  std::vector<Vertex>   vertexTable_;

  IntSet                deadEdges_;
  IntSet                deadVertices_;
};

template<typename T>
class Polygon<T>::Edge {
public:
  // Constructors
  explicit Edge() : vertex_(UNDEFINED), prevEdge_(UNDEFINED), nextEdge_(UNDEFINED), present_(true) {}
	
  Edge(index_t pVertex, index_t lpEdge, index_t lnEdge) : 
    vertex_(pVertex), prevEdge_(lpEdge), nextEdge_(lnEdge), present_(true) {}
	
  // Copy Constructor
  Edge(const Edge& he) : vertex_(he.vertex_), prevEdge_(he.prevEdge_), 
                         nextEdge_(he.nextEdge_), present_(he.present_) {}
	
  Edge& operator=(const Edge &he) {
    vertex_   = he.vertex_;
    prevEdge_ = he.prevEdge_;
    nextEdge_ = he.nextEdge_;
    present_  = he.present_;
    return *this;
  }
	
  void setVertex(index_t v) {
    vertex_ = v;
  }
	
  void setPrevEdge(index_t e) {
    prevEdge_ = e;
  }
	
  void setNextEdge(index_t e) {
    nextEdge_ = e;
  }
		
  index_t getVertex()   const {return vertex_; }
  index_t getPrevEdge() const {return prevEdge_; }
  index_t getNextEdge() const {return nextEdge_; }
	
  bool operator==(Edge& he) {
    return vertex_   == he.vertex_ 
      && prevEdge_ == he.prevEdge_ 
      && nextEdge_ == he.nextEdge_ 
      && present_  == he.present_;
  }
	
  index_t vertex_; 	
  index_t prevEdge_;  	
  index_t nextEdge_;  	
  bool present_;  	// true if the Face exists; false if deleted;
};

template<typename T>
class Polygon<T>::Vertex {
public:

  Vertex() : edge_(UNDEFINED), present_(true) {
    v_.clear();
  }
	
  Vertex(Point2d &v, int edge) : edge_(edge), present_(true) {
    v_.resize(dim);
    for (int i = 0; i < dim; i++) {
      v_[i] = static_cast<T>(v[i]);
    }
  }
	
  Vertex(PolygonVertex const &v, index_t const edge) : v_(v) , edge_(edge), present_(true) {}
	
	
  // copy constructor
  Vertex(Vertex const &v) {
    v_    = v.v_;
    edge_ = v.edge_;
    present_ = v.present_;
  };
	
  // withinEpsilon should return true if and only if the vertex is located 
  // within a distance epsilon of the position w, using an Lp metric indexed
  // by p. (The default metric is Euclidean.)	
  bool isWithinEpsilon(	PolygonVertex const &w, 
                        double const epsilon, 
                        double const p=2.0) const {
    return v_.isWithinEpsilon(w, epsilon, p);
  }
	
  bool isWithinEpsilon(	Polygon<T>::Vertex const &w, 
                        double const epsilon, 
                        double const p=2.0) const {
    return v_.isWithinEpsilon(w.v_, epsilon, p);
  }
				
  PolygonVertex	v_;		 	  // An n-tuple of type T.
  index_t		edge_;        // edgeTable_ index to the edge that follows v_;
  bool			present_;
};


// =====================
// = Private Functions =
// =====================

// TODO: improve the efficiency of the following function by searching through the
// vertex table for p, and then use the topology of the edgeTable to find n.
template<typename T>
inline bool Polygon<T>::findEdge(PolygonVertex const &p, 
                                 PolygonVertex const &n, 
                                 index_t *index, 
                                 double epsilon) const {
  for(index_t i = 0; i < edges(); i++) {
    if (vertexTable_[getCurrentVertex(i)].isWithinEpsilon(p, epsilon) && 
        edgeTable_[i].present_) {
      index_t theNextVertexIndex = getNextVertex(i);
      if (isDefined(theNextVertexIndex)) {
        if (vertexTable_[theNextVertexIndex].isWithinEpsilon(n, epsilon)) {
          *index = i;
          return true;
        }
      }
    }
  }
  return false;
}

template<typename T>
inline bool Polygon<T>::findEdge(index_t initialVertexIndex, 
                                 index_t finalVertexIndex, 
                                 index_t *index) const {
	// ensure that the first two arguments point to existing vertices.
	assert(0 <= initialVertexIndex && initialVertexIndex < vertices());
	assert(0 <= finalVertexIndex   && finalVertexIndex   < vertices());
	
    bool status = false;

	// retain the index of the Edge associated with the initial vertex.
	int edge = vertexTable_[initialVertexIndex].edge_;
    if (getNextVertex(edge) == finalVertexIndex) {
      *index = edge;
      status = true;
    }
	

	return status;  // Twin edge was not found.
}									

template<typename T>
inline bool Polygon<T>::findVertex(	PolygonVertex const &v, 
                                    index_t *index, 
                                    double epsilon) const {	
	for(index_t i = 0; i < vertices(); i++) {
		if (vertexTable_[i].isWithinEpsilon(v, epsilon)) {
			*index = i;
			return true;
		}
	}
	return false;
}


// ====================
// = Public Functions =
// ====================

template<typename T>
Polygon<T>::Polygon(std::vector<PolygonVertex> const &points) : initialEdge_(0) {
  index_t n = points.size();  // number of vertices and edges

  index_t prev = n - 1;
  for (index_t e = 0; e < n; e++) {
    index_t next = (e+1) % n;
    edgeTable_.push_back(Edge(e, prev, next));
    vertexTable_.push_back(Vertex(points[e], e));
    prev = e % n;
  }
}

template<typename T>
inline bool Polygon<T>::addEdge(PolygonVertex &v1, index_t prev, index_t next, index_t *newEdge) {
  // find vertex v1 in the vertex table
  index_t v1Index;

  if (! findVertex(v1, &v1Index, epsilonFind)) {
    // v1 was not found in vertexTable_, so add it.
    v1Index = vertexTable_.size();
    vertexTable_.push_back(Vertex(v1, UNDEFINED));
  } else if (isDefined(vertexTable_[v1Index].edge_)) {
    // cannot overwrite an existing edge
    std::cerr << "Polygon<T>::addEdge(...): "
              << "cannot overwrite an existing edge = " << vertexTable_[v1Index].edge_
              << " (referenced by vertex index = " << v1Index
              << std::endl;
    return false;
  }

  return addEdge(v1Index, prev, next, newEdge);
}

template<typename T>
bool Polygon<T>::createTriangle(PolygonVertex &v0, PolygonVertex &v1, PolygonVertex &v2) {
  if (isDefined(initialEdge_))
    return false;
  
  initialEdge_ = 0;
  vertexTable_.push_back(Vertex(v0, 0));
  vertexTable_.push_back(Vertex(v1, 1));
  vertexTable_.push_back(Vertex(v2, 2));
  edgeTable_.push_back(Edge(0, 2, 1));
  edgeTable_.push_back(Edge(1, 0, 2));
  edgeTable_.push_back(Edge(2, 1, 0));
  return true;
}

template<typename T>
inline bool Polygon<T>::addEdge(index_t vertex, index_t prevEdge, index_t nextEdge, index_t *newEdgeIndex) {
 
//  if (isDefined(prevEdge)) {
//     if (isDefined(edgeTable_[prevEdge].nextEdge_)) {
//       std::cerr << "Polygon<T>::addEdge: ERROR: " 
//                 << "prevEdge = " << prevEdge 
//                 << " has a defined nextEdge = " << edgeTable_[prevEdge].nextEdge_
//                 << std::endl;
//       return false;
//     }
//   }

//   if (isDefined(nextEdge)) {
//     if (isDefined(edgeTable_[nextEdge].prevEdge_)) {
//       std::cerr << "Polygon<T>::addEdge: ERROR: " 
//                 << "nextEdge = " << nextEdge
//                 << " has a defined prevEdge = " << edgeTable_[nextEdge].prevEdge_
//                 << std::endl;
//       return false;
//     }
//   }

  *newEdgeIndex = edgeTable_.size();
  vertexTable_[vertex].edge_ = *newEdgeIndex;
  edgeTable_.push_back(Edge(vertex, prevEdge, nextEdge));
 
  if (isDefined(prevEdge)) {
    edgeTable_[prevEdge].nextEdge_ = *newEdgeIndex;
  }

  if (isDefined(nextEdge)) {
    edgeTable_[nextEdge].prevEdge_ = *newEdgeIndex;
  }
  return true;
}

// deleteEdge marks the indicated edge for deletion and removes references to it from its
// previous and next edges. In the event that the deleted edge coincides with the initialEdge_
// of this polygon, then deleteEdges attempts to reassign initialEdge_ to an alternate
// existing edge in the same polygon.
template<typename T>
bool Polygon<T>::deleteEdge(index_t edge) {
	assert(0 <= edge && edge < edges());
	
	if (edgeTable_[edge].present_ == false) {
		std::cerr << 
			"Polygon<T>::deleteEdge(edgeIndex = " << edge << 
			"): WARNING: attempt to delete a deleted edge."	<< std::endl;
		return false;
	}
	
    int prev = getPrevEdge(edge);
    int next = getNextEdge(edge);
    int vert = getCurrentVertex(edge);

    if (isDefined(prev)) edgeTable_[prev].nextEdge_   = UNDEFINED;
    if (isDefined(next)) edgeTable_[next].prevEdge_   = UNDEFINED;
    if (isDefined(vert)) deleteVertex(vert);

	edgeTable_[edge].present_ = false;
	deadEdges_.insert(edge);
  
    if (edge == initialEdge_) {  // If true, then try to set initialEdge_ to another edge in the same polygon.
      index_t e = edge;
      do {
        e = getNextEdge(e);
        if (isPresentEdge(e)) {
          initialEdge_ = e;
        }
      } while (e != edge && initialEdge_ == edge);
      
      if (edge == initialEdge_) {
		// Attempt to reset initailEdge_ failed.
        std::cerr << "Polygon<T>::deleteEdge(" << edge << "): "
                  << "could not update initialEdge_."
                  << std::endl;
        return false;
      }
    }
	return true;
}

template<typename T>
bool Polygon<T>::deleteVertex(index_t vertexIndex) {
	assert(0 <= vertexIndex && vertexIndex < vertices());
	vertexTable_[vertexIndex].edge_ = UNDEFINED;
	vertexTable_[vertexIndex].present_ = false;
    deadVertices_.insert(vertexIndex);
	return true;
}
	


// Function purgeVertices removes from vertexTable_ all records of vertices that have
// been marked for deletion. A vertex is marked for deletion by invoking function
// deleteVertex. The latter inserts the index value of the marked vertex into the
// buffer deadVertices_.
template<typename T>
bool Polygon<T>::purgeVertices(void) {
 IntSet::iterator pos;
  for (pos = deadVertices_.begin(); pos != deadVertices_.end(); ++pos) {
    index_t lastIndex = vertexTable_.size() - 1; // the last (active) face in the table.
    if (*pos != lastIndex) {
      // Move the last face (lastIndex) into position *pos in the vertexTable_.

      // First copy the last face entries into position *pos. This purges vertex *pos.
      vertexTable_[*pos] = vertexTable_[lastIndex];

      // Then change all references to vertex lastIndex that appear in the edgeTable_
      // to vertex *pos.

      index_t edge = vertexTable_[lastIndex].edge_;
      if (isDefined(edge)) 
        edgeTable_[edge].vertex_ = *pos;
    }

    // Finally delete the last (active) entry, since it is now redundant.
    vertexTable_.pop_back();
  }

  // Once all vertices have been purged, clear the deadVertices_ buffer.
  deadVertices_.clear();
  return true;
}


template<typename T>
bool Polygon<T>::purgeEdges(void) {
 IntSet::iterator pos;
  for (pos = deadEdges_.begin(); pos != deadEdges_.end(); ++pos) {
    index_t lastIndex = edgeTable_.size() - 1;
    if (*pos != lastIndex) {
      // fill position *pos in the table with the last face.
      edgeTable_[*pos] = edgeTable_[lastIndex];
	
      // Update the vertex table, if needed.
      index_t vertexIndex = edgeTable_[*pos].vertex_;
      if (vertexTable_[vertexIndex].edge_ == lastIndex) {
		vertexTable_[vertexIndex].edge_ = *pos;
      }
	
      // Update the edge table:

      if (initialEdge_ == lastIndex) {
        initialEdge_ = *pos;
      }

      index_t nextEdge = getNextEdge(*pos);
      if (getPrevEdge(nextEdge) == lastIndex) {
		edgeTable_[nextEdge].prevEdge_ = *pos;
      } else {
		std::cerr << "Polygon<T>::deleteEdge2:: mismatch between "
                  << "edgeTable_[" << *pos << "].nextEdge_ = "
                  << edgeTable_[*pos].nextEdge_ << " and "
                  << "edgeTable_[" << nextEdge << "].prevEdge_ = "
                  << edgeTable_[nextEdge].prevEdge_;
		return false;
      }
	
      index_t prevEdge = getPrevEdge(*pos);
      if (getNextEdge(prevEdge) == lastIndex) {
		edgeTable_[prevEdge].nextEdge_ = *pos;
      } else {
		std::cerr << "Polygon<T>::deleteEdge2:: mismatch between "
                  << "edgeTable_[" << *pos << "].prevEdge_ = "
                  << edgeTable_[*pos].prevEdge_ << " and "
                  << "edgeTable_[" << prevEdge << "].nextEdge_ = "
                  << edgeTable_[prevEdge].nextEdge_;
		return false;
      }
    }
    edgeTable_.pop_back();
  }

  deadEdges_.clear();
  return true;
}


template<typename T>
bool Polygon<T>::purgeTables(void) {
 
  bool status = purgeVertices();
  if (! status) {
    std::cerr << "Polygon<T>::purgeTables: ERROR: "
              << "Failed to purgeVertices."
              << std::endl;
    return false;
  }
 
  status = purgeEdges();
  if (! status) {
    std::cerr << "Polygon<T>::purgeTables: ERROR: "
              << "Failed to purgeEdges."
              << std::endl;
    return false;
  }

  return true;
}

template<typename T>
void Polygon<T>::describe(std::ostream &os) const {
	const int formatWidth = 12;
	index_t nVertices = vertices();
	index_t nEdges    = edges();
	
    os << std::endl << std::endl;
	os << "Vertex Table (" << nVertices << ")" << std::endl;
	os 	<< std::setw(formatWidth) << "V" 
		<< std::setw(formatWidth*2) << "Position" 
		<< std::setw(formatWidth) << "Edge"
		<< std::setw(formatWidth) << "Present (T/F)" << std::endl;
		
	for(index_t i = 0; i < nVertices; ++i) {
		const Vertex &v = vertexTable_[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << std::setprecision(8) << v.v_ 
			<< std::setw(formatWidth) << v.edge_ 
			<< std::setw(formatWidth) << v.present_ << std::endl;
	}
	
	os << std::endl<< "Edge Table (" << nEdges << "), initialEdge_ = " << initialEdge_
       << std::endl;
	os 	<< std::setw(formatWidth)     << "       E" 
		<< std::setw(formatWidth+3)   << "  Vertex" 
		<< std::setw(formatWidth+1)   << "NextEdge" 
		<< std::setw(formatWidth)     << "PrevEdge" 
		<< std::setw(formatWidth+4)   << "Present (T/F)" << std::endl;
	for(index_t i = 0; i < nEdges; ++i) {
		const Edge &e = edgeTable_[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << e.vertex_ 
			<< std::setw(formatWidth) << e.nextEdge_
		   	<< std::setw(formatWidth) << e.prevEdge_ 
			<< std::setw(formatWidth) << e.present_ << std::endl;
	}
	
	os 		<< "Area      = " << this->area() << std::endl
            << "Perimeter = " << this->perimeter() << std::endl;
	
}

// isConvex(bool verbose) returns true if the polyhedron is convex; false, if it is nonconvex.
// If verbose=true, then the first nonconvex corner is described in std::cout.

template<typename T>
bool Polygon<T>::isConvex(bool verbose) const {
  index_t e1 = initialEdge_;
  if (! isDefined(e1)) abort();

  PolygonVertex v1 = vertexTable_[getCurrentVertex(e1)].v_;
 
  index_t e2 = getNextEdge(e1);
  if (! isDefined(e2)) abort();
  PolygonVertex v2    = vertexTable_[getCurrentVertex(e2)].v_;
  PolygonVertex v2mv1 = v2 - v1;

  index_t e3 = getNextEdge(e2);
  if (! isDefined(e3)) abort();

  PolygonVertex v3    = vertexTable_[getCurrentVertex(e3)].v_;
  PolygonVertex v3mv2 = v3 - v2;

    
  do {
    double prod = v2mv1 * perp(v3mv2);
    if (prod < 0) {
      if (verbose) {
        std::cout << "Polygon<T>::isConvex: First nonconvex triple: edges " 
                  << e1 << "-" << e2 << "-" << e3 
                  << ", prod = " << prod 
                  << std::endl;
      }
      return false;
    }
    v2mv1 = v3mv2;
    e1 = e2;
    e2 = e3;
    v2 = v3;
    e3 = getNextEdge(e3);
    if (! isDefined(e3)) abort();
    v3 = vertexTable_[getCurrentVertex(e3)].v_;
    v3mv2 = v3 - v2;
  } while (e1 != initialEdge_);
  return true;
}


// Function area() returns the area of the current polygon, or -1
// if the polygon is not closed.

template<typename T>
double Polygon<T>::area() const {
  double a = 0;

  index_t edge = initialEdge_;
  index_t prev = edgeTable_[edge].prevEdge_;
  do {
    a +=  vertexTable_[prev].v_ * perp(vertexTable_[edge].v_);
    prev = edge;
    edge = getNextEdge(edge);
    if (! isDefined(edge)) {
      return -1.0;
    }

  } while (edge != initialEdge_);
	
  return a/2.0;	// Inner product		
}

// Function perimeter() returns the perimeter of the current polygon, or -1
// if the polygon is not closed.
template<typename T>
double Polygon<T>::perimeter() const {
  double p = 0;

  index_t edge = initialEdge_;
  index_t prev = edgeTable_[edge].prevEdge_;
  do { 
    PolygonVertex delta = vertexTable_[edge].v_ - vertexTable_[prev].v_;
    p += delta.norm();

    prev = edge;
    edge = getNextEdge(edge);
    if (! isDefined(edge)) {
      return -1.0;
    }
  } while (edge != initialEdge_);

	return p;
}


template<typename T>
int Polygon<T>::countInactiveEdges() const {
  int count = 0;
  for(index_t i = 0; i < edgeTable_.size(); ++i) {
    if (! edgeTable_[i].present_) count++;
  }
  return count;
}

template<typename T>
int Polygon<T>::countInactiveVertices() const {
  int count = 0;
  for(index_t i = 0; i < vertexTable_.size(); ++i) {
    if (! vertexTable_[i].present_) count++;
  }
  return count;
}

template<typename T>
bool Polygon<T>::checkTables() const {
  bool status = true;

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

 
  index_t dEdges = deadEdges_.size();
  if (dEdges != xEdges) {
    std::cerr << "Warning: number of deadEdges (" << dEdges << ") "
              << "disagrees with number of inactive edges (" << xEdges << ")."
              << std::endl;
    status = false;
  }

  index_t dVertices = deadVertices_.size();
  if (dVertices != xVertices) {
    std::cerr << "Warning: number of deadVertices (" << dVertices << ") "
              << "disagrees with number of inactive vertices (" << xVertices << ")."
              << std::endl;
    status = false;
  }

  return status;
}
    
// Function topologyTest tests the topological consistency of the polygon defined
// by the tables:  edgeTable_, and vertexTable_. If the argument
// verbose is true, then a detailed error report is printed on std::cerr. If no
// topological errors are detected, then topologyTest returns true, otherwise false.

template<typename T>
bool Polygon<T>::topologyTest(bool verbose) const {
  int edgeCount = 0;
  int errorCount = 0;

  // for each active face
  for(index_t e = 0; e < edgeTable_.size(); e++) {
    if (edgeTable_[e].present_) {
      edgeCount++;

      // test next & prev consistency
      index_t nextEdge = getNextEdge(e);
      if (getPrevEdge(nextEdge) != e) {
        if (verbose) {
          std::cerr << "topologyTest: Pointer inconsistency between edge " << e
                    << " and its successor " << nextEdge 
                    << std::endl;
        }
        errorCount++;
      }
	
      // test vertex pointer consistency
      index_t vIndex = edgeTable_[e].vertex_;
      if (vertexTable_[vIndex].edge_ != e) {
        if (verbose) {
          std::cerr << "topologyTest: Pointer inconsistency between edge " << e
                    << " and its vertex " << vIndex 
                    << std::endl;
        }
        errorCount++;
      }
    }
  }
	
  index_t edge = initialEdge_;
  do {
    edgeCount--;

    edge = getNextEdge(edge);
    if (! isDefined(edge)) {
      if (verbose) {
        std::cerr << "Polygon<T>::topologyTest: encounted undefined edge during traversal."
                  << std::endl;
      }
      return false;
    }
  } while (edge != initialEdge_);

  if (edgeCount != 0) {
    errorCount++;
  }

  if (verbose && errorCount > 0) {
    std::cerr << "Polygon<T>::topologyTest:: "
              << errorCount << " errors detected."
              << std::cerr;
  }
  return (errorCount == 0);		
}
	
template<typename T>
void Polygon<T>::drawPerimeter() const {

  index_t edge = initialEdge_;
  glBegin(GL_LINE_LOOP);
  do {
    index_t v = edgeTable_[edge].vertex_;
    double x = static_cast<double>(vertexTable_[v].v_[0]);
    double y = static_cast<double>(vertexTable_[v].v_[1]);

    glVertex2d(x, y);   // OpenGL draw line segment to these coordinates.

    edge = getNextEdge(edge);
    if (! isDefined(edge)) {
      std::cerr << "Polygon<T>::drawPerimeter(): "
                << " polygon is incomplete."
                << std::endl;
    }
  } while (edge != initialEdge_);
  glEnd();
}
  
#endif // __POLYGON_H__
