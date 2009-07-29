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
#include <iomanip>
#include <valarray>

const int dim = 3;

typedef double Point3d[3];
typedef float  Point3f[3];
typedef long   Point3i[3];

// A polyhedron is an ordered set of three-dimensional vertices, edges, and 
// faces. The typename T designates the type of each vertex component (e.g.,
// char, short, int, float, or double.)
template<typename T>
class Polyhedron {
public:

  // typedef double normalType;
  typedef ntuple<double,3>  FaceNormal;
  typedef ntuple<T,3>       FaceVertex;
  typedef std::set<int, std::greater<int> >  IntSet;

  static const int SEARCH_FAILED = -1;
	
  Polyhedron() {} // default constructor
	
  // copy constructor
  Polyhedron(Polyhedron<T> const &p) {
    edgeTable_   = p.edgeTable_;
    vertexTable_ = p.vertexTable_;
    faceTable_   = p.faceTable_;
  }
	
  // addFace adds a new face, defined by the given vertex list, to the 
  // current polyhedral model. The index of the created face is returned.
  bool addFace(std::list<FaceVertex> const &vList, int* faceIndex);
	
  bool addFace(FaceVertex const &a, FaceVertex const &b, FaceVertex const &c, int* faceIndex) {
    std::list<FaceVertex> face;
    face.push_back(a);
    face.push_back(b);
    face.push_back(c);
    return addFace(face, faceIndex);	
  }
	
  // addFace(int, int, int) adds a triangular face with the indicated
  // indexed vertices. The index of the created face is returned.
  bool addFace(int v1Index, int v2Index, int v3Index, int* faceIndex);
	
  // marks the indicated face (along with its implicit edges, and vertices,
  // if they are unique) for deletion from the current model.
  // deleteFace returns true if the face was deleted without detecting any
  // errors. If false is returned, then an error was detected, and the 
  // calling routine should either throw an exception or abort.
  bool deleteFace(int faceIndex);
  bool deleteFace2(int faceIndex);
	
  // marks the indicated edge for deletion.
  bool deleteEdge(int edgeIndex);
  bool deleteEdge2(int edgeIndex);
	
  // marks the indicated vertex for deletion.
  bool deleteVertex(int vertexIndex);
  bool deleteVertex2(int vertexIndex);
	
  // Extends the face that contains the adjEdgeIndex by adding edges between 
  // the vertex with index v1Index and the endpoints of the edge adjEdgeIndex.
  // If successful, the boolean value true is returned. Otherwise, false
  // is returned.
  bool extendFace(int v1Index, int adjEdgeIndex);
	
  // getCollocalEdge searches for an edge different from edgeIndex that shares its vertex.
  // If found, the index of the discovered collocal edge is returned via the pointer in
  // the second argument, and the function returns the value true. If no edge can be
  // found, then false is returned. Note that three search strategies are employed in
  // succession.
  bool getCollocalEdge(int edgeIndex, int *collocalEdge) const;
	
  // removes all faces, halfedges, and vertices that are marked for deletion.
  // purgeFaces returns the number of faces that were purged.
  bool purgeFaces();
  bool purgeEdges();
  bool purgeVertices();
  bool purgeTables();
		
	
  // Returns the number of faces in the current polytope.
  int faces() const    {return faceTable_.size();  }
	
  // Returns the number of edges in the current polytope.
  int edges() const    {return edgeTable_.size();  }
	
  // Returns the number of vertices in the current polytope.
  int vertices() const {return vertexTable_.size();}
	
  // isConvex() returns true of the polyhedron is convex.
  bool isConvex() const;
	
  // Returns the next available index for the vertex table.
  int nextVertexIndex() {return vertices();}
	
  // Returns the next available index for the edge table.
  int nextEdgeIndex() {return edges();}
	
  // Returns the next available index for the edge table.
  int lastEdgeIndex() {return edges() - 1;}
	
  // Returns the next available index for the face table.
  int nextFaceIndex() {return faces();}
	
  // Returns the number of edges that lack a twinEdge.
  int boundaryEdges() const {
    int count = 0;
    for(int i = 0; i < edgeTable_.size(); ++i) {
      if (edgeTable_[i].twinEdge_ == UNDEFINED) 
        ++count;
    }
    return count;
  }
	
  FaceNormal &getNormalVector(int faceIndex) {
    assert(0 <= faceIndex && faceIndex < faces());
    return faceTable_[faceIndex].normalVector_;
  }
	
  FaceVertex &getFaceVertex(int faceIndex) {
    assert(0 <= faceIndex && faceIndex < faces());
    return vertexTable_[edgeTable_[faceTable_[faceIndex].outer_].vertex_].v_;
  }

  // A vertex v is said to be exterior to face f, if v lies on the
  // positive side of the plane that coincides with f. Vertex v, is said to be
  // coplanar with respect to f, if v lies in the plane that coincides with f,
  // and is said to be interior to face f, if v lies on the negative side
  // of the plane that coincides with face f. The function visibilityTest(...) 
  // returns +1 if the vertex v is exterior, 0 if it is coplanar, and -1 if v
  // is interior to f. For queries of floating point data, the third argument
  // epsilon can be used to define a slab of with epsilon that is centered
  // on face f.
  int visibilityTest(int faceIndex, FaceVertex& v, T epsilon=0) {
  }

  // Returns the number of edges that have a twinEdge.
  int interiorEdges() const {
    return edges() - boundaryEdges();
  }
	
  // Checks the three tables (faceTable_, edgeTable_, and vertexTable_) for the presence of
  // inactive elements. If the first argument is set to true, then a detailed error report
  // is printed on std::cerr. If no inactive elements are detected, then checkTables returns
  // true, otherwise false.
  bool checkTables(bool verbose) const;

  // Function countInactiveFaces returns the number of faces in the faceTable_ that are not
  // active, i.e. have present_ set to false.
  int countInactiveFaces() const;

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
  bool eulerTest(bool verbose = false) const; 
	
  // Prints a description of the current polyhedron on the indicated stream.
  void describe(std::ostream &os) const;
	
  // Computes are returns the area of the face with index i.
  double faceArea(int i) const;
 
  // Computes and returns the surface are of the current polyhedron.
  double surfaceArea() const ;

  // Computes and returns the volume of the current polyhedron.
  double volume() const;
	
protected:
  class HalfEdge;  // Nested class for representing each HalfEdge
  class Vertex;    // Nested class for representing each Vertex
  class Face;      // Nested class for representing each Face
	
  // Serves as a null pointer.
  static const int UNDEFINED = -1;
  static const int poly_level = 1; // debug level
	
#if (DEBUG 	> 0)	
  static const bool debug = (DEBUG >= poly_level);
#else
  static const bool debug = false;
#endif
	
  int getUndefined() const {return UNDEFINED;}
	
  // Returns true if a halfedge already exists that extends from vertex p 
  // to vertex n, within precision epsilon. If true, then index returns the
  // location of this edge in the HalfEdge table.
  bool findEdge(FaceVertex const &p, FaceVertex const &n, int *index, 
                double epsilon = 0) const;
	
  // findEdge(int, int, int) returns true if there is an edge in the edge
  // table that extends from the initial vertex to the final vertex. In
  // this case, the index of the discovered edge is returned via the 
  // third argument. If no edge is found, false is returned.
  bool findEdge(int initialVertexIndex, int finalVertexIndex, int *index) const;
	
  // Returns true if a vertex v already exists in the vertex table, within 
  // precision epsilon. If true, then index returns the location of this 
  // vertex in the vertex table.
  bool findVertex(FaceVertex const &v, int* index, double epsilon = 0) const;
	
  // Returns the index of the current Face for the indicated edge.
  int getCurrentFace(int edge) const {
    return (isDefined(edge) ? edgeTable_[edge].face_ : getUndefined());
  }
	
  // Returns the index of the current Vertex for the indicated edge.
  int getCurrentVertex(int edge) const {
    return (isDefined(edge) ? edgeTable_[edge].vertex_ : getUndefined());
  }
	
  // Returns the halfedge index that follows the value edge.
  int getNextEdge(int edge) const {
    return (isDefined(edge) ? edgeTable_[edge].nextEdge_ : getUndefined());
  }
	
	
  int getClockwiseEdge(int edge) const {
    return getNextEdge(getTwinEdge(edge));
  }
	
  int getCounterClockwiseEdge(int edge) const {
    return getTwinEdge(getPrevEdge(edge));
  }
	
  // Returns the index of the next Vertex for the indicated edge.
  int getNextVertex(int edge) const {
    return getCurrentVertex(getNextEdge(edge));
  }
	
  // Returns the previous half edge
  int getPrevEdge(int edge) const {
    return (isDefined(edge) ? edgeTable_[edge].prevEdge_ : getUndefined());
  }
		
  // Returns the halfedge index that precedes the value edge.
  int getPrevVertex(int edge) const {
    return getCurrentVertex(getPrevEdge(edge));
  }
	
  // Returns the halfedge index that is opposite to the value edge.
  int getTwinEdge(int edge) const {
    return (isDefined(edge) ? edgeTable_[edge].twinEdge_ : getUndefined()); 
  }
	
  // Returns true if the halfedge with index value edge has a twin.
  bool hasTwin(int edge) const {
    return (getTwinEdge(edge) != UNDEFINED);
  }
	
  bool isDefined(int index) const {
    return (index != UNDEFINED);
  }

  // Function isolatedVertexCheckAndDelete deletes the edge with value edgeToDelete after
  // checking to see if its vertex is used by any other edge in the edge circuit that 
  // contains the edge with index successor. If no other edge has the same vertex, then
  // the vertex is also deleted. Otherwise, the edge field of the vertex entry is updated
  // if necessary, to ensure that there exists a reference to an active edge. Also, the
  // outer edge field of the face entry is updated, if necessary. Note that edgeToDelete
  // and successor should belong (formerly) to the same face.
  bool isolatedVertexCheckAndDelete(int edgeToDelete, int successor);

  bool isPresentEdge(int index) const {
    return (0 <= index && index < edges() && edgeTable_[index].present_);
  }
	
  int edgeOfVertex(int vertexIndex) {
    return (isDefined(vertexIndex) ? vertexTable_[vertexIndex].edge_ : getUndefined());
  }
	
  bool updateEdgesOfVertex(int vIndex, int vIndexOld);
	
  bool vertexIsIsolated(int vertexIndex) {
    return edgeOfVertex(vertexIndex) == UNDEFINED; 
  }
	
  // Compute the normal vector to the indicated face.
  void computeNormalVector(int faceIndex);
	
  // Compute the normal vectors of every face.
  void computeNormalVectors();
	
  // Merges the two indicated faces if they are adjacent and coplanar.
  bool mergeFaces(int fold, int fnew);
	
  // Clear the polytope.
  void clear() {
    edgeTable_.clear();
    vertexTable_.clear();
    faceTable_.clear();
  }
	
  // ================
  // = Data members =
  // ================
	
  std::vector<HalfEdge> edgeTable_;
  std::vector<Vertex>   vertexTable_;
  std::vector<Face>	    faceTable_;

  IntSet    deadFaces_;
  IntSet    deadEdges_;
  IntSet    deadVertices_;
};

template<typename T>
class Polyhedron<T>::HalfEdge {
public:
  // Constructors
  explicit HalfEdge() : vertex_(UNDEFINED), prevEdge_(UNDEFINED), 
                        nextEdge_(UNDEFINED), twinEdge_(UNDEFINED), face_(UNDEFINED), present_(true) {}
	
  HalfEdge(int pVertex, int lpEdge, int lnEdge, int twinEdge, int lFace) : 
    vertex_(pVertex), prevEdge_(lpEdge), nextEdge_(lnEdge), 
    twinEdge_(twinEdge), face_(lFace), present_(true) {}
	
  // Copy Constructor
  HalfEdge(const HalfEdge& he) : vertex_(he.vertex_), 
                                 prevEdge_(he.prevEdge_), nextEdge_(he.nextEdge_), twinEdge_(he.twinEdge_),
                                 face_(he.face_), present_(he.present_) {}
	
  HalfEdge& operator=(const HalfEdge &he) {
    vertex_   = he.vertex_;
    prevEdge_ = he.prevEdge_;
    nextEdge_ = he.nextEdge_;
    twinEdge_ = he.twinEdge_;
    face_     = he.face_;
    present_  = he.present_;
    return *this;
  }
	
  void setVertex(int pVertex) {
    vertex_ = pVertex;
  }
	
  void setPrevEdge(int lpEdge) {
    prevEdge_ = lpEdge;
  }
	
  void setNextEdge(int lnEdge) {
    nextEdge_ = lnEdge;
  }
	
  void setTwinEdge(int twinEdge) {
    twinEdge_ = twinEdge;
  }
	
  void setFace(int lFace) {
    face_ = lFace;
  }
	
  int getVertex()   const {return vertex_; }
  int getPrevEdge() const {return prevEdge_; }
  int getNextEdge() const {return nextEdge_; }
  int getTwinEdge() const {return twinEdge_;}
  int getFace()     const {return face_; }
	
  bool operator==(HalfEdge& he) {
    return vertex_   == he.vertex_ 
      && prevEdge_ == he.prevEdge_ 
      && nextEdge_ == he.nextEdge_ 
      && twinEdge_ == he.twinEdge_ 
      && face_     == he.face_
      && present_  == he.present_;
  }
	
  int vertex_; 	
  int prevEdge_;  	
  int nextEdge_;  
  int twinEdge_; 	
  int face_;	
  bool present_;  	// true if the Face exists; 
  // false if deleted;
};

template<typename T>
class Polyhedron<T>::Face {
public:

	
  Face() : outer_(UNDEFINED), present_(true) {
    normalVector_.clear();
  }
	
  explicit Face(int outerEdge) : outer_(outerEdge), present_(true) {
    normalVector_.clear();
  }
	
  void setOuterEdge(int e) {outer_ = e;}
	
  int 				outer_;
  std::vector<int> 	inner_;	
  FaceNormal          normalVector_;
  bool 				present_;  	// true if the Face exists; false if deleted;
};


template<typename T>
class Polyhedron<T>::Vertex {
public:

  Vertex() : edge_(UNDEFINED), present_(true) {
    v_.clear();
  }
	
  Vertex(Point3d &v, int edge) : edge_(edge), present_(true) {
    v_.resize(dim);
    for (int i = 0; i < dim; i++) {
      v_[i] = static_cast<T>(v[i]);
    }
  }
	
  Vertex(FaceVertex const &v, int const edge) : v_(v) , edge_(edge), present_(true) {}
	
	
  // copy constructor
  Vertex(Vertex const &v) {
    v_    = v.v_;
    edge_ = v.edge_;
    present_ = v.present_;
  };
	
  // withinEpsilon should return true if and only if the vertex is located 
  // within a distance epsilon of the position w, using an Lp metric indexed
  // by p. (The default metric is Euclidean.)	
  bool isWithinEpsilon(	FaceVertex const &w, 
                        double const epsilon, 
                        double const p=2.0) const {
    return v_.isWithinEpsilon(w, epsilon, p);
  }
	
  bool isWithinEpsilon(	Polyhedron<T>::Vertex const &w, 
                        double const epsilon, 
                        double const p=2.0) const {
    return v_.isWithinEpsilon(w.v_, epsilon, p);
  }
				
  FaceVertex 	v_;		  // An n-tuple of type T.
  int 			edge_;    // HalfEdge* edge_;
  bool			present_;
};


// =====================
// = Private Functions =
// =====================

// TODO: improve the efficiency of the following function by searching through the
// vertex table for p, and then use the topology of the edgeTable to find n.
template<typename T>
inline bool Polyhedron<T>::findEdge(FaceVertex const &p, 
									FaceVertex const &n, 
									int *index, 
									double epsilon) const {
  for(int i = 0; i < edges(); i++) {
    if (vertexTable_[getCurrentVertex(i)].isWithinEpsilon(p, epsilon) && 
        edgeTable_[i].present_) {
      int theNextVertexIndex = getNextVertex(i);
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
inline bool Polyhedron<T>::findEdge(	int initialVertexIndex, 
										int finalVertexIndex, 
										int *index) const {
	// ensure that the first two arguments point to existing vertices.
	assert(0 <= initialVertexIndex && initialVertexIndex < vertices());
	assert(0 <= finalVertexIndex   && finalVertexIndex   < vertices());
	
	// retain the index of the HalfEdge associated with the initial vertex.
	int initialEdgeIndex = vertexTable_[initialVertexIndex].edge_;
	
	if (! isPresentEdge(initialEdgeIndex)) {
		// std::cerr << "Polyhedron<T>::findEdge: initial edge for vertex " << initialVertexIndex
		// 			<< " is not present in the edge table." << std::endl;
		return false;
	}
	
	// As we walk around the vertex (clockwise), edge indicates the current
	// edge that points away from the initial vertex.
	int nEdges = edges();
	int edge = initialEdgeIndex;
	do {
		if (getNextVertex(edge) == finalVertexIndex) {  // Eureka!
			*index = edge; // return the discovered edge index.
			return true;
		}
		
		// Search the edges that eminanate from the initial vertex in a clockwise order.
		edge = getNextEdge(getTwinEdge(edge));
	} while(edge != initialEdgeIndex && isPresentEdge(edge));
	
	if (edge != initialEdgeIndex) {	
		// Starting again from the initialEdge for the initial vertex, search the edges that
		// eminate from the initial vertex in a counter-clockwise order.
		
		initialEdgeIndex = getTwinEdge(getPrevEdge(edge));  // since the initial edge was checked
															// previously, skip forward.
		if (isPresentEdge(initialEdgeIndex)) {
			int edge = initialEdgeIndex;
			do {
				if (getNextVertex(edge) == finalVertexIndex) {  // Eureka!
					*index = edge; // return the discovered edge index.
					return true;
				}

				edge = getTwinEdge(getPrevEdge(edge)); // search counter-clockwise about initial vertex.
			
				// Continue walking until either the twin is undefined, or we
				// return to the initalEdgeIndex.
			} while(edge != initialEdgeIndex && isPresentEdge(edge));	
		}
	}
	
	// Now search around the final vertex
	initialEdgeIndex = getPrevEdge(vertexTable_[finalVertexIndex].edge_);
	if (! isPresentEdge(initialEdgeIndex)) {
		initialEdgeIndex = getTwinEdge(vertexTable_[finalVertexIndex].edge_);
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
			int edge = getPrevEdge(getTwinEdge(initialEdgeIndex)); // counter-clockwise step
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
										int *index, 
										double epsilon) const {	
	for(int i = 0; i < vertices(); i++) {
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
inline bool Polyhedron<T>::addFace(typename std::list< FaceVertex > const &v, int* faceIndex) {
	typename std::list< FaceVertex >::const_iterator preVertex;
	
	int size = v.size(); 				// number of new edges to add to edgeTable_
	int edgeIndexBase   = edges();		// index of next new HalfEdge	
	int edgeIndex       = edgeIndexBase;
	*faceIndex          = faces();		// index of next new Face
	
	Face newFace(edgeIndex);
	faceTable_.push_back(newFace);
			
	for(preVertex = v.begin(); preVertex != v.end(); ++preVertex, ++edgeIndex) {
		int vertexIndex;		// index of next new Vertex
		
		int nextEdgeIndex = edgeIndex + 1;
		if (nextEdgeIndex == edgeIndexBase + size) // the nextEdgeIndex of the untimate edge is edgeIndexBase.
			nextEdgeIndex = edgeIndexBase;
		
		int prevEdgeIndex = edgeIndex - 1;
		if (prevEdgeIndex == edgeIndexBase - 1) // the prevEdgeIndex of the first edge is (edgeIndexBase + size - 1).
			prevEdgeIndex += size;
		
		if (! findVertex(*preVertex, &vertexIndex)) {
			// vertex *preVertex is not in the vertexTable_, therefore, ...
			vertexIndex = vertices();                   // set vertexIndex to the next available vertex index.
			Vertex newVertex(*preVertex, edgeIndex);	// create a new vertex for *preVertex,
			vertexTable_.push_back(newVertex);          // add the new vertex to vertexTable_.
		}
		
		// Find the twin edge, if it exists. First step: find the next vertex
		// in the list of points.
		typename std::list< FaceVertex >::const_iterator nextVertex = preVertex;
		if (distance(preVertex, v.end()) > 1) {
			++nextVertex; 
		} else {
			// for the last edge, the next vertex is the initial one.
			nextVertex = v.begin();
		}
		
		// The twin edge should point from the next vertex to the current one.
		int twinEdgeIndex;
		if (findEdge(*nextVertex, *preVertex, &twinEdgeIndex)) {
			if (edgeTable_[twinEdgeIndex].twinEdge_ != UNDEFINED) {
				std::cerr << 
					"Polyhedron<T>::addFace: ERROR: current face may conincide with a previous face."
					<< std::endl;
				return false;
			}
			edgeTable_[twinEdgeIndex].twinEdge_ = edgeIndex;
		} else {
			twinEdgeIndex = UNDEFINED;
		}
		
		HalfEdge newEdge(vertexIndex, prevEdgeIndex, nextEdgeIndex, 
						twinEdgeIndex, *faceIndex);
		edgeTable_.push_back(newEdge);
	}
	
	computeNormalVector(*faceIndex);
	
	return true;
}

// Function isolatedVertexCheckAndDelete is an auxilliary function used by function enlargeFace
// (see below).The first argument, edgeIndex, should identify a half-edge (E) that belongs to a
// closed polygonal face (F) contained in the internal tables. If the vertex (V) of this half-edge
// does not belong to any other half-edge in F, then both the V and E will be deleted by inserting
// them in the deadVertices_ and deadEdges_ sets respectively. If V belongs to another half-edge,
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
bool Polyhedron<T>::isolatedVertexCheckAndDelete(int edgeToDelete, int goodEdge) {
  int vIndex = getCurrentVertex(edgeToDelete);
  int altEdge = UNDEFINED;
  int edge = goodEdge;

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
      if (vertexTable_[vIndex].edge_ == edgeToDelete) {
        vertexTable_[vIndex].edge_ = altEdge;
      }
      break;
    }
          
    edge = getNextEdge(edge);
  } while (edge != goodEdge);

  if (! isDefined(altEdge)) {   // No edge was found
    deadVertices_.insert(vIndex);
    vertexTable_[vIndex].present_ = false;
  }

  int faceIndex = getCurrentFace(edgeToDelete);
  if (faceTable_[faceIndex].outer_ == edgeToDelete) {
    faceTable_[faceIndex].outer_ = goodEdge;
  }

  deadEdges_.insert(edgeToDelete);
  edgeTable_[edgeToDelete].present_ = false;
  return true;
}


// Function extendFace adds a coplanar triangular region (R) to an existing face (F). Two
// arguments are provided that define the R and F: v1Index is the index of vertex V1 that 
// belongs to R, and adjEdgeIndex is the index of an edge, directed from vertex V3 to V2,
// that belongs to F (but does not contian V1). Thus R is defined as the triangle that 
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
inline bool Polyhedron<T>::extendFace(int v1Index, int adjEdgeIndex) {
	assert(0 <= v1Index && v1Index < vertices());               // Validate arguments.
	assert(0 <= adjEdgeIndex && adjEdgeIndex < edges());
	
	int twin21;                                                 // v2->v1
	int twin13;                                                 // v1->v3
	int v2Index   = getNextVertex(adjEdgeIndex);                // v2
	int v3Index   = getCurrentVertex(adjEdgeIndex);             // v1
	int faceIndex = getCurrentFace(adjEdgeIndex);               // F
	

	if (findEdge(v2Index, v1Index, &twin21)) {                  // Does v2->v1 exist?
		int oldTwin = edgeTable_[twin21].twinEdge_;             // If so, it should be twinless.
		if (isDefined(oldTwin) && edgeTable_[oldTwin].present_) {
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
		int oldTwin = edgeTable_[twin13].twinEdge_;             // If so, it should be twinless.
		if (isDefined(oldTwin) && edgeTable_[oldTwin].present_) {
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
		int successor = getNextEdge(twin21);
    
        if (! isolatedVertexCheckAndDelete(twin21, successor)) {
          std::cerr << "Polyhedron<T>::extendFace(" << v1Index << ", "
                    << adjEdgeIndex << "): " 
                    << "Undefined edge encounted while traversing face "
                    << faceIndex << std::endl;
          return false;
        }

		if (getCurrentFace(twin13) == faceIndex) {              // Case 4?
			// The extension fills in a triangular hole.
			int antecedent = getPrevEdge(twin13);
			edgeTable_[antecedent].nextEdge_ = successor;
			edgeTable_[successor ].prevEdge_ = antecedent;
              
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
			edgeTable_[adjEdgeIndex].nextEdge_ = successor;
			edgeTable_[adjEdgeIndex].twinEdge_ = twin13;
			if (isDefined(twin13)) {
				edgeTable_[twin13].twinEdge_ = adjEdgeIndex;
			}
			edgeTable_[successor].prevEdge_ = adjEdgeIndex;
		}
	} else if (getCurrentFace(twin13) == faceIndex) {           // Case 3?
		// The extension fills a notch between vertices v2, v3, v1.
		int successor = getNextEdge(adjEdgeIndex);
		edgeTable_[twin13].nextEdge_ = successor;
		edgeTable_[successor].prevEdge_ = twin13;
		edgeTable_[twin13].twinEdge_ = twin21;
		if (isDefined(twin21)) {
			edgeTable_[twin21].twinEdge_ = twin13;
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
		int nextEdgeIndex = edges();
		int successor = getNextEdge(adjEdgeIndex);
		edgeTable_[adjEdgeIndex].nextEdge_ = nextEdgeIndex;
		edgeTable_[adjEdgeIndex].twinEdge_ = twin13;
		edgeTable_[successor].prevEdge_ = nextEdgeIndex;
		if (isDefined(twin13)) {
			edgeTable_[twin13].twinEdge_ = adjEdgeIndex;
		}
		if (isDefined(twin21)) {
			edgeTable_[twin21].twinEdge_ = nextEdgeIndex;
		}
		
		HalfEdge edge12(v1Index, adjEdgeIndex, successor, twin21, faceIndex);
		edgeTable_.push_back(edge12);
		if (! isDefined(vertexTable_[v1Index].edge_)) 
			vertexTable_[v1Index].edge_ = nextEdgeIndex;
	}
	
	return true;
}


template<typename T>
inline bool Polyhedron<T>::addFace(int v1Index, int v2Index, int v3Index, int* faceIndex) {
	assert(0 <= v1Index && v1Index < vertices());
	assert(0 <= v2Index && v2Index < vertices());
	assert(0 <= v3Index && v3Index < vertices());
	
	int twin21 = UNDEFINED;
	int twin13 = UNDEFINED;
	int twin32 = UNDEFINED;
		
	int nextEdgeIndex = edges();
	*faceIndex = faces();
	
	if (findEdge(v2Index, v1Index, &twin21)) {
		int oldTwin = edgeTable_[twin21].twinEdge_;
		if (isDefined(oldTwin) && edgeTable_[oldTwin].present_) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin21
						<< " connecting vertex " << v2Index << " to " << v1Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
		//edgeTable_[twin21].twinEdge_ = nextEdgeIndex;
	} else {
		twin21 = UNDEFINED;
	}
	
	if (findEdge(v3Index, v2Index, &twin32)) {
		int oldTwin = edgeTable_[twin32].twinEdge_;
		if (isDefined(oldTwin) && edgeTable_[oldTwin].present_) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin32
						<< " connecting vertex " << v3Index << " to " << v2Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
		//edgeTable_[twin32].twinEdge_ = nextEdgeIndex + 1;
	} else {
		twin32 = UNDEFINED;
	}
	
	if (findEdge(v1Index, v3Index, &twin13)) {
		int oldTwin = edgeTable_[twin13].twinEdge_;
		if (isDefined(oldTwin) && edgeTable_[oldTwin].present_) {
			std::cerr 	<< "Polyhedron<T>::addFace(" << v1Index << ", "
						<< v2Index << ", " << v3Index << ", ...): Found twin edge " << twin13
						<< " connecting vertex " << v1Index << " to " << v3Index
						<< " that is bound to twin " << oldTwin  
						<< std::endl;
			return false;
		}
		// edgeTable_[twin13].twinEdge_ = nextEdgeIndex + 2;
	} else {
		twin13 = UNDEFINED;
	}
	
	if (isDefined(twin21)) edgeTable_[twin21].twinEdge_ = nextEdgeIndex;
	if (isDefined(twin32)) edgeTable_[twin32].twinEdge_ = nextEdgeIndex + 1;
	if (isDefined(twin13)) edgeTable_[twin13].twinEdge_ = nextEdgeIndex + 2;
	
	if (vertexIsIsolated(v1Index)) {
		vertexTable_[v1Index].edge_ = nextEdgeIndex;
	}
	
	if (vertexIsIsolated(v2Index)) {
		vertexTable_[v2Index].edge_ = nextEdgeIndex + 1;
	}
	
	if (vertexIsIsolated(v3Index)) {
		vertexTable_[v3Index].edge_ = nextEdgeIndex + 2;
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
					
	edgeTable_.push_back(edge12);
	edgeTable_.push_back(edge23);
	edgeTable_.push_back(edge31);
	
	Face newFace(nextEdgeIndex);
	faceTable_.push_back(newFace);
	computeNormalVector(*faceIndex);
	return true;
}

template<typename T>
void Polyhedron<T>::computeNormalVector(int faceIndex) {
	Face* face = &faceTable_[faceIndex];
	FaceNormal &normal = face->normalVector_;
	
	int edge = face->outer_;
	int next = edgeTable_[edge].nextEdge_;
	normal.clear();
	
	do {
		FaceVertex& p = vertexTable_[getCurrentVertex(edge)].v_;
		FaceVertex& n = vertexTable_[getNextVertex(edge)].v_;
					
		normal[0] += (p[2] + n[2])*(p[1] - n[1]);
		normal[1] += (p[0] + n[0])*(p[2] - n[2]);
		normal[2] += (p[1] + n[1])*(p[0] - n[0]);
		
		edge = next;
		next = edgeTable_[edge].nextEdge_;
	} while (edge != face->outer_);
	
	// normalize the normal vector
	normal.normalize();
}

template<typename T>
void Polyhedron<T>::computeNormalVectors() {
	for(int faceIndex = 0; faceIndex < faces(); ++faceIndex) {
		computeNormalVector(faceIndex);
	}
}


// deleteFace marks the indicated face for removal, and functionally deletes
// the face from the table, including its edges, and any vertices that exclusively
// belong to deleted faces. 

template<typename T>
bool Polyhedron<T>::deleteFace2(int faceIndex) {
	assert(0 <= faceIndex && faceIndex < faces());
	
	// Identify the edges that delineate the selected face, and store them in the
	// list edgeLoop.
	int initialEdge = faceTable_[faceIndex].outer_;
	int edge = initialEdge;

    // The sets deadEdges and deadVertices are used to accumulate the indices
    // of edges and vertices to be deleted in decreasing order. (The default,
    // which uses std::less<int>, constructs the sets in increasing order.)
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

        int vertexIndex = getCurrentVertex(edge);
		if (isDefined(vertexIndex)) {
			if (vertexTable_[vertexIndex].edge_ == edge) {
				int altEdge;
				bool status = getCollocalEdge(edge, &altEdge);
				if (status) {
					vertexTable_[vertexIndex].edge_ = altEdge;
				} else {
					deadVertices.insert(vertexIndex);
				}
			}
		}

		edge = getNextEdge(edge);	
	} while(edge != initialEdge);
	
	assert(edge == initialEdge);
	do {
		int vertexIndex = getCurrentVertex(edge);
		if (isDefined(vertexIndex)) {
			if (vertexTable_[vertexIndex].edge_ == edge) {
				int altEdge;
				bool status = getCollocalEdge(edge, &altEdge);
				if (status) {
					vertexTable_[vertexIndex].edge_ = altEdge;
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
	
	int lastFace = faces() - 1;
	faceTable_[faceIndex] = faceTable_[lastFace];
	faceTable_.pop_back();
	
	initialEdge = faceTable_[faceIndex].outer_;
	edge = initialEdge;
	do {
		if (! isDefined(edge)) {
			std::cerr << "Warning[deleteFace2(" << faceIndex << ")]: "
					  << "relocated an incomplete face."
					  << std::endl;
			return false;
		}
		edgeTable_[edge].face_ = faceIndex;
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
bool Polyhedron<T>::getCollocalEdge(int edgeIndex, int *collocalEdge) const {
	int vertexIndex = edgeTable_[edgeIndex].vertex_;
	*collocalEdge = getUndefined();
	
	if (isDefined(vertexIndex)) {
		*collocalEdge = getClockwiseEdge(edgeIndex); // returns an edge that is collocal with edgeIndex
												     // or UNDEFINED.
		if (! isDefined(*collocalEdge) || (! edgeTable_[*collocalEdge].present_)) {
			*collocalEdge = getCounterClockwiseEdge(edgeIndex);
		}
		
		if (! isDefined(*collocalEdge) || (! edgeTable_[*collocalEdge].present_)) {
			for(int edge = 0; edge < edges(); ++edge) {
				if (edgeTable_[edge].vertex_ == vertexIndex && 
					edge != edgeIndex &&
					edgeTable_[edge].present_) {
						*collocalEdge = edge;
						break;
					}
			}	
		}
	}
	return isDefined(*collocalEdge) && edgeTable_[*collocalEdge].present_;
}

// delete the indicated edge and remove it from the edgeTable_ by overwriting
// edgeTable_[edgeIndex] with the last item in the table, and then popping the
// edgeTable_. Before calling this routine, all references to the edge to be
// removed should have been removed from the faceTable_ and vertexTable_.
// Since the index of the last edge in the table has been changed from edges()
// to edgeIndex, all references to the last edge in the faceTable_, edgeTable_,
// and vertexTable_ should be updated.

template<typename T>
bool Polyhedron<T>::deleteEdge2(int edgeIndex) {
	// Delete the indicated edge by overwriting
	int lastEdgeIndex = edges() - 1;
	edgeTable_[edgeIndex] = edgeTable_[lastEdgeIndex];
	edgeTable_.pop_back();
	
	// Update the face table, if needed.
	int faceIndex = edgeTable_[edgeIndex].face_;
	if (faceTable_[faceIndex].outer_ == lastEdgeIndex) {
		faceTable_[faceIndex].outer_ = edgeIndex;
	}
	
	// Update the vertex table, if needed.
	int vertexIndex = edgeTable_[edgeIndex].vertex_;
	if (vertexTable_[vertexIndex].edge_ == lastEdgeIndex) {
		vertexTable_[vertexIndex].edge_ = edgeIndex;
	}
	
	// Update the edge table:

	int nextEdge = getNextEdge(edgeIndex);
	if (getPrevEdge(nextEdge) == lastEdgeIndex) {
		edgeTable_[nextEdge].prevEdge_ = edgeIndex;
	} else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
			<< "edgeTable_[" << edgeIndex << "].nextEdge_ = "
			<< edgeTable_[edgeIndex].nextEdge_ << " and "
			<< "edgeTable_[" << nextEdge << "].prevEdge_ = "
			<< edgeTable_[nextEdge].prevEdge_;
		return false;
	}
	
	int twinEdge = getTwinEdge(edgeIndex);
	if (getTwinEdge(twinEdge) == lastEdgeIndex) {
		edgeTable_[twinEdge].twinEdge_ = edgeIndex;
	} else if (isDefined(getTwinEdge(twinEdge))) {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
			<< "edgeTable_[" << edgeIndex << "].twinEdge_ = "
			<< edgeTable_[edgeIndex].twinEdge_ << " and "
			<< "edgeTable_[" << twinEdge << "].twinEdge_ = "
			<< edgeTable_[twinEdge].twinEdge_;
		return false;		
	}
	
	int prevEdge = getPrevEdge(edgeIndex);
	if (getNextEdge(prevEdge) == lastEdgeIndex) {
		edgeTable_[prevEdge].nextEdge_ = edgeIndex;
	} else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
			<< "edgeTable_[" << edgeIndex << "].prevEdge_ = "
			<< edgeTable_[edgeIndex].prevEdge_ << " and "
			<< "edgeTable_[" << prevEdge << "].nextEdge_ = "
			<< edgeTable_[prevEdge].nextEdge_;
		return false;
	}
	
	return true;
}

// deleteVertex2 deletes the indicated vertex from vertexTable_ by replacing 
// entry vIndex with the last entry in the table. Before this procedure is
// called, all edges that reference vertex vIndex should have been deleted.
// All edges that radiate from the relocated vertex (those that are collocal),
// need to be updated with the new index.
template<typename T>
bool Polyhedron<T>::deleteVertex2(int vIndex) {
	int lastVertex = vertices() - 1;
	vertexTable_[vIndex] = vertexTable_[lastVertex];
	vertexTable_.pop_back();
	
	return Polyhedron<T>::updateEdgesOfVertex(vIndex, lastVertex);	
}

// Function updateEdgesOfVertex replaces each reference to vertex vOld to the new value
// vNew in the entries in edgeTable_. It is assumed that all edges that eminate from vertex
// vOld are topologically accessible by pointer algebra.

template<typename T>
bool Polyhedron<T>::updateEdgesOfVertex(int vNew, int vOld) {
  int initialEdge = vertexTable_[vOld].edge_;
  assert(0 <= vNew < vertices());
  assert(0 <= vOld < vertices());


  if (isDefined(initialEdge)) {		
    enum State {CLOCKWISE, COUNTERCLOCKWISE, SEQUENTIAL, FINISHED};
    State state = CLOCKWISE;
    int edge = initialEdge;
    
    while (state != FINISHED) {
      switch (state) {
      case CLOCKWISE:
        if (isPresentEdge(edge)) {
          edgeTable_[edge].vertex_ = vNew;
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
          edgeTable_[edge].vertex_ = vNew;
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
          if (edgeTable_[edge].vertex_ == vOld && isPresentEdge(edge)) {
            edgeTable_[edge].vertex_ = vNew;
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
bool Polyhedron<T>::deleteFace(int faceIndex) {
	assert(0 <= faceIndex && faceIndex < faces());
	
	if (faceTable_[faceIndex].present_ == false) {
		std::cerr << 
			"Polyhedron<T>::deleteFace(faceIndex = " << faceIndex << 
			"): WARNING: attempt to delete a deleted face."	<< std::endl;
		return false;
	}
	
	faceTable_[faceIndex].present_ = false;  // marks the face for deletion
	deadFaces_.insert(faceIndex);

	int initialEdge = faceTable_[faceIndex].outer_;
	if (isDefined(initialEdge)) {
		int edge = initialEdge;
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
bool Polyhedron<T>::deleteEdge(int edgeIndex) {
	assert(0 <= edgeIndex && edgeIndex < edges());
	
	if (edgeTable_[edgeIndex].present_ == false) {
		std::cerr << 
			"Polyhedron<T>::deleteEdge(edgeIndex = " << edgeIndex << 
			"): WARNING: attempt to delete a deleted edge."	<< std::endl;
		return false;
	}
	
	edgeTable_[edgeIndex].present_ = false;
	deadEdges_.insert(edgeIndex);

	// remove any reference to edge by its twin.
	int twinEdge = getTwinEdge(edgeIndex);
	if (isDefined(twinEdge)) {
		if (edgeTable_[twinEdge].twinEdge_ == edgeIndex) {
			edgeTable_[twinEdge].twinEdge_ = UNDEFINED; 
		} else {
			std::cerr <<
				"Polyhedron<T>::deleteEdge(edgeIndex = " << edgeIndex << 
				"): ERROR: twin inconsistency with edge "
				<< twinEdge << std::endl;
            return false;
		}
	}		
	
    // Modify faceTable_ references to edgeIndex:
	int faceIndex = edgeTable_[edgeIndex].face_;
	if (faceTable_[faceIndex].outer_ == edgeIndex && 
		faceTable_[faceIndex].present_) {
		// If possible, replace faceTable_[faceIndex].outer_ with
		// an edge that is retained (or present).
		
		int edge = getNextEdge(edgeIndex);
		do {
			if (edgeTable_[edge].face_ == faceIndex &&
				edgeTable_[edge].present_) {
					faceTable_[faceIndex].outer_ = edge;
					break;
			}
			edge = getNextEdge(edge);
		} while (edge != edgeIndex && isDefined(edge));

        if (faceTable_[faceIndex].outer_ == edgeIndex) {
          // Could not update faceTable_[faceIndex].outer_ with an active edge.
          faceTable_[faceIndex].outer_ = UNDEFINED;
        }
	}
	
    // Modify vertexTable_ references to edgeIndex:
	int vertexIndex = edgeTable_[edgeIndex].vertex_;
	if (vertexTable_[vertexIndex].edge_ == edgeIndex) {
		// If possible, replace vertexTable_[edgeTable_[edge].vertex_].edge_
		// with an active edge that originates from vertexIndex.
		int edge = getClockwiseEdge(edgeIndex);
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
					} else if (edgeTable_[edge].present_) {
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
					} else if (edgeTable_[edge].present_) {
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
						if (edgeTable_[edge].vertex_ == vertexIndex && edgeTable_[edge].present_) {
							state = SUCCESS;
							break;
						}
					}
					if (state != SUCCESS) {
						state = FAILURE;
					}
					break;
					
				case SUCCESS:
					vertexTable_[vertexIndex].edge_ = edge;
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
bool Polyhedron<T>::deleteVertex(int vertexIndex) {
	assert(0 <= vertexIndex && vertexIndex < vertices());
	vertexTable_[vertexIndex].edge_ = UNDEFINED;
	vertexTable_[vertexIndex].present_ = false;
    deadVertices_.insert(vertexIndex);
	return true;
}
	
// purgeFaces removes all faces that are marked for deletion, i.e., those
// faces that contain a value of present_ == false. In addition, all halfEdges
// and vertices that are dependent on a face that is removed are deleted.
template<typename T>
bool Polyhedron<T>::purgeFaces(void) {
  IntSet::iterator pos;
  for (pos = deadFaces_.begin(); pos != deadFaces_.end(); ++pos) {
    int lastFaceIndex = faceTable_.size() - 1;
    if (*pos != lastFaceIndex) {
      // fill position *pos in the table with the last face.
      faceTable_[*pos] = faceTable_[lastFaceIndex];
      // update cross references from edgeTable
      int initialEdge = faceTable_[*pos].outer_;
      int edge = initialEdge;
      do {
		if (! isDefined(edge)) {
          std::cerr << "Polyhedron<T>::purgeFaces: "
                    << "faceTable_[" << lastFaceIndex << "] is incomplete."
                    << std::endl;
          return false;
		}

		edgeTable_[edge].face_ = *pos;
		edge = getNextEdge(edge);
      } while (edge != initialEdge);
    }
    faceTable_.pop_back();
  }
  deadFaces_.clear();
  return true;
}

// Function purgeVertices removes from vertexTable_ all records of vertices that have
// been marked for deletion. A vertex is marked for deletion by invoking function
// deleteVertex. The latter inserts the index value of the marked vertex into the
// buffer deadVertices_.
template<typename T>
bool Polyhedron<T>::purgeVertices(void) {
 IntSet::iterator pos;
  for (pos = deadVertices_.begin(); pos != deadVertices_.end(); ++pos) {
    int lastIndex = vertexTable_.size() - 1; // the last (active) face in the table.
    if (*pos != lastIndex) {
      // Move the last face (lastIndex) into position *pos in the vertexTable_.

      // First copy the last face entries into position *pos. This purges vertex *pos.
      vertexTable_[*pos] = vertexTable_[lastIndex];

      // Then change all references to vertex lastIndex that appear in the edgeTable_
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
    vertexTable_.pop_back();
  }

  // Once all vertices have been purged, clear the deadVertices_ buffer.
  deadVertices_.clear();
  return true;
}


template<typename T>
bool Polyhedron<T>::purgeEdges(void) {
 IntSet::iterator pos;
  for (pos = deadEdges_.begin(); pos != deadEdges_.end(); ++pos) {
    int lastIndex = edgeTable_.size() - 1;
    if (*pos != lastIndex) {
      // fill position *pos in the table with the last face.
      edgeTable_[*pos] = edgeTable_[lastIndex];

      // Update the face table, if needed.
      int faceIndex = edgeTable_[*pos].face_;
      if (faceTable_[faceIndex].outer_ == lastIndex) {
		faceTable_[faceIndex].outer_ = *pos;
      }
	
      // Update the vertex table, if needed.
      int vertexIndex = edgeTable_[*pos].vertex_;
      if (vertexTable_[vertexIndex].edge_ == lastIndex) {
		vertexTable_[vertexIndex].edge_ = *pos;
      }
	
      // Update the edge table:

      int nextEdge = getNextEdge(*pos);
      if (getPrevEdge(nextEdge) == lastIndex) {
		edgeTable_[nextEdge].prevEdge_ = *pos;
      } else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
                  << "edgeTable_[" << *pos << "].nextEdge_ = "
                  << edgeTable_[*pos].nextEdge_ << " and "
                  << "edgeTable_[" << nextEdge << "].prevEdge_ = "
                  << edgeTable_[nextEdge].prevEdge_;
		return false;
      }
	
      int twinEdge     = getTwinEdge(*pos);
      int twinTwinEdge = getTwinEdge(twinEdge);
      if (isDefined(twinTwinEdge)) {
        if (twinTwinEdge == lastIndex) {
          edgeTable_[twinEdge].twinEdge_ = *pos;
        } else {
          std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
                    << "edgeTable_[" << *pos << "].twinEdge_ = "
                    << edgeTable_[*pos].twinEdge_ << " and "
                    << "edgeTable_[" << twinEdge << "].twinEdge_ = "
                    << edgeTable_[twinEdge].twinEdge_;
          return false;		
        }
      }
	
      int prevEdge = getPrevEdge(*pos);
      if (getNextEdge(prevEdge) == lastIndex) {
		edgeTable_[prevEdge].nextEdge_ = *pos;
      } else {
		std::cerr << "Polyhedron<T>::deleteEdge2:: mismatch between "
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
	int nVertices = vertices();
	int nEdges    = edges();
	int nFaces    = faces();
	
	os << "Vertex Table (" << nVertices << ")" << std::endl;
	os 	<< std::setw(formatWidth) << "V" 
		<< std::setw(formatWidth*2) << "Position" 
		<< std::setw(formatWidth) << "Edge"
		<< std::setw(formatWidth) << "Present (T/F)" << std::endl;
		
	for(int i = 0; i < nVertices; ++i) {
		const Vertex &v = vertexTable_[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << std::setprecision(8) << v.v_ 
			<< std::setw(formatWidth) << v.edge_ 
			<< std::setw(formatWidth) << v.present_ << std::endl;
	}
	
	os << std::endl<< "Edge Table (" << nEdges << ")" << std::endl;
	os 	<< std::setw(formatWidth)     << "       E" 
		<< std::setw(formatWidth+3)   << "  Vertex" 
		<< std::setw(formatWidth+1)   << "NextEdge" 
		<< std::setw(formatWidth)     << "PrevEdge" 
		<< std::setw(formatWidth)     << "TwinEdge" 
		<< std::setw(formatWidth-3)   << "    Face" 
		<< std::setw(formatWidth+4)   << "Present (T/F)" << std::endl;
	for(int i = 0; i < nEdges; ++i) {
		const HalfEdge &e = edgeTable_[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << e.vertex_ 
			<< std::setw(formatWidth) << e.nextEdge_
		   	<< std::setw(formatWidth) << e.prevEdge_ 
			<< std::setw(formatWidth) << e.twinEdge_ 
			<< std::setw(formatWidth) << e.face_
			<< std::setw(formatWidth) << e.present_ << std::endl;
	}
	
	os 	<< std::endl << "Face Table (" << nFaces << ")" << std::endl;
	os 	<< std::setw(formatWidth) << "F" 
		<< std::setw(formatWidth+4) << "OuterEdge" 
		<< std::setw(formatWidth+4)   << "NormalVector"
		<< std::setw(formatWidth+4)   << "Present (T/F)" << std::endl;
	
	for(int i = 0; i < nFaces; ++i) {
		const Face &f = faceTable_[i];
		
		os 	<< std::setw(formatWidth) << i 
			<< std::setw(formatWidth) << f.outer_  
			<< std::setw(formatWidth) << f.normalVector_ 
			<< std::setw(formatWidth) << f.present_<< std::endl;
	}
	
	os 	<< "Volume = " << this->volume() 
		<< "; Surface Area = " << this->surfaceArea() << std::endl;
	
}

// isConvex() returns true if the polyhedron is convex; false, if it is nonconvex.

template<typename T>
bool Polyhedron<T>::isConvex() const {

	typename std::vector<typename Polyhedron<T>::Face>::const_iterator facePos;
	for(facePos = faceTable_.begin(); facePos != faceTable_.end(); ++facePos) {
		if (facePos->present_) {
			FaceVertex faceCenter; // The center of the current face.
			int vertexCount;		// Number of vertices that define the current face.

			// Compute the face center.
			faceCenter.clear();  // Initialize accumulators.
			vertexCount = 0;

			int initialEdge = facePos->outer_;  // Identify the initial edge.
			int edge = initialEdge;
			do {	// for each edge that delimits the current face.
				faceCenter += vertexTable_[getCurrentVertex(edge)].v_;
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
			for(vertexPos = vertexTable_.begin(); vertexPos != vertexTable_.end(); ++vertexPos) {
				FaceVertex displacement = vertexPos->v_ - faceCenter;
				double ip = facePos->normalVector_ * displacement.normalize();  // inner product
				if (ip  > 0) {
					std::cout << "Nonconvex relation: " 
						<< facePos->normalVector_ << "* ("
						<< vertexPos->v_ << " - " 
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
double Polyhedron<T>::faceArea(int faceIndex) const {
	FaceVertex edgeProd;
	int initialEdge = faceTable_[faceIndex].outer_;
	int edge = initialEdge;
	do {
		edgeProd += crossProduct(vertexTable_[getCurrentVertex(edge)].v_, 
								 vertexTable_[getNextVertex(edge)].v_);
		edge = getNextEdge(edge);
	} while (edge != initialEdge);
	
	return (edgeProd*faceTable_[faceIndex].normalVector_)/2.0;	// Inner product		
}

template<typename T>
double Polyhedron<T>::surfaceArea() const {
	T area = 0;
	for (int i = 0; i < faces(); i++) {
		if (faceTable_[i].present_) {
			area += faceArea(i);
		}
	}
	return area;
}

template<typename T>
int Polyhedron<T>::countInactiveFaces() const {
  int count = 0;
  for(int i = 0; i < faceTable_.size(); ++i) {
    if (! faceTable_[i].present_) count++;
  }
  return count;
}

template<typename T>
int Polyhedron<T>::countInactiveEdges() const {
  int count = 0;
  for(int i = 0; i < edgeTable_.size(); ++i) {
    if (! edgeTable_[i].present_) count++;
  }
  return count;
}

template<typename T>
int Polyhedron<T>::countInactiveVertices() const {
  int count = 0;
  for(int i = 0; i < vertexTable_.size(); ++i) {
    if (! vertexTable_[i].present_) count++;
  }
  return count;
}

template<typename T>
bool Polyhedron<T>::checkTables(bool verbose) const {
  bool status = true;

  int xFaces = countInactiveFaces();
  if (xFaces > 0) {
    std::cerr << "Warning: " << xFaces << " are inactive." << std::endl;
    status = false;
  }

  int xEdges = countInactiveEdges();
  if (xEdges > 0) {
    std::cerr << "Warning: " << xEdges << " are inactive." << std::endl;
    status = false;
  }

  int xVertices = countInactiveVertices();
  if (xVertices > 0) {
    std::cerr << "Warning: " << xVertices << " are inactive." << std::endl;
    status = false;
  }

  int dFaces = deadFaces_.size();
  if (dFaces != xFaces) {
    std::cerr << "Warning: number of deadFaces (" << dFaces << ") "
              << "disagrees with number of inactive faces (" << xFaces << ")."
              << std::endl;
    status = false;
  }

  int dEdges = deadEdges_.size();
  if (dEdges != xEdges) {
    std::cerr << "Warning: number of deadEdges (" << dEdges << ") "
              << "disagrees with number of inactive edges (" << xEdges << ")."
              << std::endl;
    status = false;
  }

  int dVertices = deadVertices_.size();
  if (dVertices != xVertices) {
    std::cerr << "Warning: number of deadVertices (" << dVertices << ") "
              << "disagrees with number of inactive vertices (" << xVertices << ")."
              << std::endl;
    status = false;
  }

  return status;
}
    
// Function eulerTest tests the topological consistency of the polyhedron defined
// by the three tables: faceTable_, edgeTable_, and vertexTable_. If the argument
// verbose is true, then a detailed error report is printed on std::cerr. If no
// topological errors are detected, then eulerTest returns true, otherwise false.

template<typename T>
bool Polyhedron<T>::eulerTest(bool verbose) const {
  std::valarray<int> eFreq(0, edgeTable_.size());
  std::valarray<int> vFreq(0, vertexTable_.size());
  int faceCount = 0;
  int edgeCount = 0;
  int vertexCount = 0;
  int errorCount = 0;
	
  // for each active face
  for(int i = 0; i < faceTable_.size(); i++) {
    if (faceTable_[i].present_) {
      faceCount++;
			
      int firstEdge = faceTable_[i].outer_;
      int edge      = firstEdge;
      int twin;
      do {

        if (! isDefined(edge) ) {
          std::cerr << "eulerTest: ERROR : face " << i << " is incomplete!"
                    << " Cannot continue test!"
                    << std::endl;
          return false;
        }

        if (! edgeTable_[edge].present_) {
          if (verbose) {
            std::cerr << "eulerTest: Edge " << edge << ", required by face " 
                      << i << " is marked for deletion." 
                      << std::endl;
          }
          errorCount++;
        }

        if ( edgeTable_[edge].face_ != i) {
          if (verbose) {
            std::cerr <<  "eulerTest: Edge " << edge << ", belongs to face " 
                      << i << " but references face " << edgeTable_[edge].face_
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
                      << " has edge " << edgeTable_[twin].twinEdge_ 
                      << " as its twin!" << std::endl;
          errorCount++;
        }

        // test next & prev consistency
        int nextEdge = getNextEdge(edge);
        if (getPrevEdge(nextEdge) != edge) {
          if (verbose) {
            std::cerr << "eulerTest: Pointer inconsistency between edge " << edge
                      << " and its successor " << nextEdge 
                      << std::endl;
          }
          errorCount++;
        }

				
        eFreq[edge]++;                                  // indicate that edge was used.
        vFreq[edgeTable_[edge].vertex_]++;              // indicate that vertex was used
				
    		
        edge = nextEdge;
      } while (edge != firstEdge);
    }
  }
	
  // Report on edge usage.
  for(int i = 0; i < eFreq.size(); i++) {
    if (eFreq[i] == 1) {
      edgeCount++;
    } else if (eFreq[i] > 1) {
      if (verbose) {
        std::cerr << "eulerTest: WARNING: edge " << i << "was used " 
                  << eFreq[i] << " times."
                  << std::endl;
      }
      errorCount++;
    } else if (edgeTable_[i].present_) {
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
	
  for(int i = 0; i < vFreq.size(); i++) {
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
  int i = 0;
	
  FaceVertex centroid;
  int vertexCount = 0;
	
  typename std::vector<Vertex>::const_iterator vertexPos;
  for(vertexPos = vertexTable_.begin(); vertexPos != vertexTable_.end(); ++vertexPos) {
    centroid += vertexPos->v_;
    vertexCount++;
  }
  centroid /= vertexCount;
	
  typename std::vector<Face>::const_iterator facePos;
  for(facePos = faceTable_.begin(); facePos != faceTable_.end(); ++facePos) {
    if (facePos->present_) {
      FaceVertex edgeProd;
      FaceVertex faceCenter;
      int count = 0;
		
      edgeProd.clear();
      faceCenter.clear();
		
      int initialEdge = facePos->outer_;
      int edge = initialEdge;
      do {
        count++;
        faceCenter += vertexTable_[getCurrentVertex(edge)].v_;
        edgeProd   += crossProduct(	vertexTable_[getCurrentVertex(edge)].v_, 
                                    vertexTable_[getNextVertex(edge)].v_);
        edge = getNextEdge(edge);
      } while (edge != initialEdge);
			
      faceCenter /= count;
			
      double sectionVolume = edgeProd*(faceCenter - centroid); // inner product	
      volume += sectionVolume;
    }
  }
  return volume/6.0; // Geometric factor (1/2)*(1/3)
}
	
#endif // __POLYHEDRON_H__
