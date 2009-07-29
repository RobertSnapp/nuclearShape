#ifndef __CONVEX_HULL_H__
#define __CONVEX_HULL_H__

//
//  convexHull.h
//
//  Created by Robert R. Snapp on 2007-02-18.
//  Copyright © 2007 Robert R. Snapp. All rights reserved.
//

#include <algorithm>  // for swap(a, b)
#include <cassert>
#include <cmath>
#include <cstdlib>  // for rand()
#include <vector>
#include <deque>
#include <list>
#include <set>
#include "polyhedron.h"
#include "ntuple.h"

const double epsilonPerp = 1.0E-08;
const double initialConditionFactor = 0.10;
const double minimumConditionFactor = 1.0E-06;

const int r_debug_max = 10;
bool const debug = false;
inline bool ch_debug(int r) {
  return false; //r==14696; //(r<=9); // 957 || r==14696);
}

template <typename T>
class ConvexHull : public Polyhedron<T> {
public:
  ConvexHull();
  ConvexHull(std::vector<typename Polyhedron<T>::FaceVertex> &points);
  bool presortConvexVector(std::vector<typename Polyhedron<T>::FaceVertex> &points, 
							 double conditionFactor = 0);
  void compute(std::vector<typename Polyhedron<T>::FaceVertex> &points);
protected:
  // static const double epsilonPerp = 1.0E-8;  
  // two unit vectors are classified as orthogonal if
  // the magnitude of their inner product is less than
  // epsilonPerp.

  void insertFacePointConflict(index_t f, index_t p);
  bool isConflictGraphConsistent(index_t npts) const;
  bool isaFacePointConflict(index_t f, index_t p) const;
  bool addFace(index_t v1, index_t v2, index_t v3, index_t* faceIndex);
  bool purgeConflictGraph(void);
  void displayConflictGraph(std::ostream &os, size_t maxIndex) const;

  std::vector<std::set<index_t> > d_faceConflicts;
  std::vector<std::set<index_t> > d_pointConflicts;
};

template<typename T>
void ConvexHull<T>::insertFacePointConflict(index_t f, index_t p) {
  d_faceConflicts[p].insert(f);
  d_pointConflicts[f].insert(p);
}

template<typename T>
bool ConvexHull<T>::isConflictGraphConsistent(index_t npts) const {
  if (Polyhedron<T>::faces() != d_pointConflicts.size()) {
    std::cerr << "ConvexHull<T>::isConflictGraphConsistent: "
              << "faces() = " << Polyhedron<T>::faces() 
              << "but d_pointConflicts.size() = " << d_pointConflicts.size()
              << std::endl;
    return false;
  }

  if (npts != d_faceConflicts.size()) {
    std::cerr << "ConvexHull<T>::isConflictGraphConsistent: "
              << "npts = " << npts 
              << "but d_faceConflicts.size() = " << d_faceConflicts.size()
              << std::endl;
    return false;
  }

  typename std::set<index_t>::iterator p;
  for (index_t i = 0; i < Polyhedron<T>::faces(); i++) {
    for(p = d_pointConflicts[i].begin(); p != d_pointConflicts[i].end(); ++p) {
      if (d_faceConflicts[*p].find(i) == d_faceConflicts[*p].end()) {
        std::cerr << "ConvexHull<T>::isConflictGraphConsistent: "
                  << "edge inconsistency between face i = " << i
                  << " and point *p = " << *p
                  << std::endl;
        return false;
      }
    }
  }

  typename std::set<index_t>::iterator f;
  for (index_t i = 0; i < npts; i++) {
    for(f = d_faceConflicts[i].begin(); f != d_faceConflicts[i].end(); ++f) {
      if (d_pointConflicts[*f].find(i) == d_pointConflicts[*f].end()) {
        std::cerr << "ConvexHull<T>::isConflictGraphConsistent: "
                  << "edge inconsistency between point i = " << i
                  << " and face *p = " << *f
                  << std::endl;
        return false;
      }
    }
  }

  return true;
}
  
template<typename T>
bool ConvexHull<T>::isaFacePointConflict(index_t f, index_t p) const {
  return (d_faceConflicts[p].find(f) != d_faceConflicts[p].end());
}


template<typename T>
bool ConvexHull<T>::addFace(index_t v1, index_t v2, index_t v3, index_t* faceIndex) {
  bool status;
  d_pointConflicts.push_back(std::set<index_t>());

  status = Polyhedron<T>::addFace(v1, v2, v3, faceIndex);
  if (Polyhedron<T>::faces() != d_pointConflicts.size()) {
    std::cerr << "ConvexHull<T>::addFace(" << v1 << ", " << v2 << ", " << v3 << ",...): "
              << "d_pointConflicts.size() = " << d_pointConflicts.size()
              << " but faces() = " << Polyhedron<T>::faces() 
              << std::endl;
    return false;
  }

  return status;
}

template<typename T>
bool ConvexHull<T>::purgeConflictGraph(void) {
  int status = true;

  typename Polyhedron<T>::IntSet::iterator f;
  for (f = Polyhedron<T>::d_deadFaces.begin(); f != Polyhedron<T>::d_deadFaces.end(); ++f) {
    index_t fLast = d_pointConflicts.size() - 1;

    // Remove references to *f from the faceConflict_ link sets.
    typename std::set<index_t>::iterator p;
    for(p = d_pointConflicts[*f].begin(); p != d_pointConflicts[*f].end(); ++p) {
      d_faceConflicts[*p].erase(*f);
    }

    if (*f != fLast) {
      d_pointConflicts[*f] = d_pointConflicts[fLast];
     
      for(p = d_pointConflicts[fLast].begin(); p != d_pointConflicts[fLast].end(); ++p) {
        if (d_faceConflicts[*p].erase(fLast)) {
          d_faceConflicts[*p].insert(*f);
        } else {
          status = false;
        }
      }
    }
    d_pointConflicts.pop_back();
  }

  return status;
}

template<typename T>
void ConvexHull<T>::displayConflictGraph(std::ostream &os, size_t maxIndex) const {
  os << "Conflict Graph: (" 
     << d_pointConflicts.size() << " face nodes, "
     << d_faceConflicts.size()  << " point nodes)"
     << std::endl;
  
  for(index_t i = 0; i < std::min(d_pointConflicts.size(), maxIndex); i++) {
    os << "L" << i << ": ";
    index_t j = 0;
    std::set<index_t>::iterator p;
    for(p = d_pointConflicts[i].begin(); p != d_pointConflicts[i].end() && j < maxIndex + 1; ++p, ++j) {
      os << *p << ", ";
      if (j == maxIndex) {
        os << "...";
      } 
    }
    os << std::endl;
  }

  os << "-----------------"
     << std::endl;

  for(index_t i = 0; i < std::min(d_faceConflicts.size(), maxIndex); i++) {
    os << "R" << i << ": ";
    index_t j = 0;
    std::set<index_t>::iterator f;
    for(f = d_faceConflicts[i].begin(); f != d_faceConflicts[i].end() && j < maxIndex + 1; ++f, ++j) {
      os << *f << ", ";
      if (j == maxIndex) {
        os << "...";
      } 
    }
    os << std::endl;
  }
}

// presortConvexVector takes a vector of three-dimensional vertices, and 
// reorders them such that the first four define a parallelpiped with positive
// volume V that satisfies the condition that V exceed the product of 
// conditionFactor times the products of its three sides. Note that 
// conditionFactor should be less than one.
//
// If successful, presortConvexVector returns true; if a parallelpiped with
// non-zero volume cannot be found, as is the case if all vectors are coplanar,
// then false is returned.
template <typename T>
bool ConvexHull<T>::presortConvexVector(
                         std::vector<typename Polyhedron<T>::FaceVertex> &points, 
                         double conditionFactor) {
	assert(0 <= conditionFactor && conditionFactor < 1.0);
	index_t size = points.size();
	
	// Randomize the vector points:
	for(index_t count = 0; count < size; count++) {
		// subtracting and adding 4 ensures that p and q will
		// skip over the first four elements.
		index_t p = (std::rand() % size);
		index_t q = (std::rand() % size);
		std::swap(points[p], points[q]);
	}
	
	for(index_t i0 = 0; i0 < size - 3; ++i0)
		for(index_t i1 = i0 + 1; i1 < size - 2; ++i1)
			for(index_t i2 = i1 + 1; i2 < size - 1; ++i2)
				for(index_t i3 = i2 + 1; i3 < size; ++i3) {
					typename Polyhedron<T>::FaceVertex v1 = points[i1] - points[i0]; // sides of the 
					typename Polyhedron<T>::FaceVertex v2 = points[i2] - points[i0]; // parallelpiped
					typename Polyhedron<T>::FaceVertex v3 = points[i3] - points[i0];
					
					T bound = conditionFactor*v1.norm()*v2.norm()*v3.norm();
					if (Polyhedron<T>::debug) {
						std::cerr << "bound = " << bound;
					}
					
					T det   = tripleProduct(v1, v2, v3); // oriented volume
					if (Polyhedron<T>::debug) {
						std::cerr << ", tripleProduct(v1, v2, v3) = " << det << std::endl;
					}
					
					if (std::abs(det) > bound) {
						if (Polyhedron<T>::debug) {
							std::cerr << "i0 = " << i0 << 
						            ", i1 = " << i1 <<
									", i2 = " << i2 <<
									", i3 = " << i3 << std::endl;
						}
						std::swap(points[0], points[i0]);
						std::swap(points[1], points[i1]);
						std::swap(points[2], points[i2]);
						std::swap(points[3], points[i3]);
						
						if (det < 0) {
							std::swap(points[1], points[2]);
						}
						
						v1 = points[1] - points[0]; // sides of the parallepiped
						v2 = points[2] - points[0];
						v3 = points[3] - points[0];
						
						assert(tripleProduct(v1, v2, v3) > 0);
						return true; // Exits with success
					}
				}
	return false; // Exits with failure.
}

template <typename T>
ConvexHull<T>::ConvexHull() : Polyhedron<T>::Polyhedron() {}

template <typename T>
ConvexHull<T>::ConvexHull(std::vector<typename Polyhedron<T>::FaceVertex > &points) {
  compute(points);
}

template <typename T>
void ConvexHull<T>::compute(std::vector<typename Polyhedron<T>::FaceVertex > &points) {
  this->clear();
  bool isConsistent = true;  // This flag is used for debugging. It should be set to false
							   // if a topological inconsistency is detected by eulerTest.
  
  index_t coplanarAdditions = 0; // counts the number of vertices that when added to the convex hull, are coplanar
                                 // with an existing face.

  index_t noncoplanarAdditions = 0; // counts the number of vertices that when added to the convex hull, are not
                                    // coplanar with an existing face.
	

  // Ensure that the model is clear.
  // Polyhedron<T>::clear();
	
  //// The following data structures have been converted into class data members.
  // The i-th element of faceConflicts should return the list of face-indicies
  // that conflict with point[i].
  //std::vector<std::set<index_t> > faceConflicts;
	
  // The j-th element of pointConflicts should return the list of point-indices
  // that conflict with face[j]
  //std::vector<std::set<index_t> > pointConflicts;
	
		
  if (points.size() < 4) {
	std::cerr << "ConvexHull<T>::ConvexHull(std::vector<FaceVertex>): "
			  << "Number of input points must exceed 3." << std::endl;
	exit(1);
  }
	
  // presortConvexVector reorders the vector points so that the first four elements form a
  // non-degenerate tetrahedron with a condition factor not less than that of the second
  // argument (conditionFactor).

  double conditionFactor = initialConditionFactor;
  bool status;
  do {
	status = presortConvexVector(points, conditionFactor);
	conditionFactor /= 2.0;
  } while (! status && conditionFactor > minimumConditionFactor);

  if (! status) {
	std::cerr << "Could not find an initial tetrahedron with conditionFactor greater than " 
			  << minimumConditionFactor << ": " << std::endl
			  << "points may be coplanar." << std::endl;
	exit(1);
  }
	
  // Construct the initial tetrahedron and verify that edge orientations are correct.
  // the first face is a triangle delemited by points[0], points[1], points[3]
	
  index_t f; // face index returned by addFace	
  status = Polyhedron<T>::addFace(points[0], points[1], points[3], &f);
  assert(status && Polyhedron<T>::getNormalVector(f)*(points[2]-points[1]) < 0);
	
  // the second face is a triangle delemited by points[1], points[2], points[3].
  status = Polyhedron<T>::addFace(points[1], points[2], points[3], &f);
  assert(status && Polyhedron<T>::getNormalVector(f)*(points[0]-points[1]) < 0);
	
  // the third face is a triangle delmited by points[0], points[3], points[2].
  status = Polyhedron<T>::addFace(points[0], points[3], points[2], &f);
  assert(status && Polyhedron<T>::getNormalVector(f)*(points[1]-points[0]) < 0);
	
  // the fourth face is a triangle delmited by points[0], points[2], points[1]
  status = Polyhedron<T>::addFace(points[0], points[2], points[1], &f);
  assert(status && Polyhedron<T>::getNormalVector(f)*(points[3]-points[0]) < 0);
	
  assert(Polyhedron<T>::faces()    ==  4);
  assert(Polyhedron<T>::edges()    == 12);
  assert(Polyhedron<T>::vertices() ==  4);
	
  if (debug) this->describe(std::cout);
	
  // The conflict graph is a bipartite graph that consists of two types of
  // nodes. On set of nodes, pointConflicts, contains one node per polyhedron
  // face. Thus, pointConflicts[i] refers to the set of vertices that 
  // conflict with face i. The other set of nodes, faceConflicts, contains
  // one node per ntuple in points. Thus, faceConflicts[j] refers to the 
  // set of faces that conflict with points[j].
	
  // Initialize the conflict graph:
  if (debug) {
	std::cout << "Creating the conflict graph ...";
  }

  // Create an initially empty set of face conflicts for each vertex contained in points.
  d_faceConflicts.resize(points.size(), std::set<index_t>());

  // Create an intially empty set for point conflicts for each of the current four (4)  faces
  // in the polyhedral model.
  d_pointConflicts.resize(4, std::set<index_t>());
   

  index_t conflicts = 0;
  for(index_t f = 0; f < 4; ++f) {
	typename Polyhedron<T>::FaceNormal faceNormal = Polyhedron<T>::getNormalVector(f);
	typename Polyhedron<T>::FaceVertex faceVertex = Polyhedron<T>::getFaceVertex(f);	
	
	
	for(index_t p = 4; p < points.size(); ++p) {
			if(faceNormal*(points[p] - faceVertex) > 0) {
				// Face f and p conflict (both cannot belong to the CH).
				conflicts++;
				insertFacePointConflict(f,p);
			}
		}
	}
		
	if (debug) {
      std::cout << "done: " <<  conflicts << " edges.\n";
      if (isConflictGraphConsistent(points.size())) {
			std::cout << "face-point conflict graph passes initial symmetry test." << std::endl;
      } else {
			std::cout << "face-point conflict graph fails initial symmetry test." << std::endl;
      }
	}
	
	// Incrementally, construct the convex hull (CH).
	for(index_t r = 4; r < points.size(); ++r) {
		// remove all faces that are in conflict with points[r]
		
		// Logarithmic "poor man's" progress bar.
		if (r %    10 == 0 && r < 250) std::cout << ".";
		if (r %   100 == 0 && r < 1000) std::cout << ","; 
		if (r %  1000 == 0 && r < 10000 && r >= 1000) std::cout << ":";
		if (r % 10000 == 0 && r >= 10000) std::cout << "!";
		
		if (d_faceConflicts[r].size()  > 0) {
			// point[r] is outside the current convex polygon and must
			// be added to the convex hull.
			if (debug && ch_debug(r)) {
				std::cout 	<< "Adding points[" << r << "]:" << std::endl
							<< "points[" << r << "] = " << points[r] << std::endl;
				std::cout   << "coplanar faces added = " << coplanarAdditions
					<< ", noncoplanar faces added = " << noncoplanarAdditions << std::endl << std::endl;
				Polyhedron<T>::describe(std::cout);
				std::cout << std::endl;
			}
		

            // conflictingFaces contains the set of face indices that conflict with points[r].
     
            std::set<index_t> conFaces = d_faceConflicts[r]; 
            std::set<index_t>::iterator pos; 

			// Identify the outermost edge loop that delemits the set of conflicting
			// faces. Note that because the polygon is by construction always convex,
			// this outerLoop delimits the portion of the polygon that is visible
			// from point[r].
			std::list<index_t> outerLoop;
			
			if (debug && ch_debug(r)) {
				std::cout << "Constructing the outer loop for points[" << r << "] ...\n";
			}
			
			// The initial edge of the outer loop is found by tracing the perimeters of each
			// face that conflicts with points[r] until a half edge is found that is adjacent
			// to (that is, the twinEdge of) a face that does not conflict with points[r].
			
			index_t initialEdge = Polyhedron<T>::getUndefined();
			
			for(pos = conFaces.begin(); 
                pos != conFaces.end() && (! Polyhedron<T>::isDefined(initialEdge)); // stop once an initial edge has been found.
				++pos) {	
				
				index_t initialFaceEdge = Polyhedron<T>::d_faceTable[*pos].d_outer;
				index_t trialEdge = initialFaceEdge;
														
				if (debug && ch_debug(r)) {
					std::cout << "trialEdge = " << trialEdge << std::endl;
				}
				
				do {
					index_t twinFace = Polyhedron<T>::getCurrentFace(
											Polyhedron<T>::getTwinEdge(trialEdge));
					if (isaFacePointConflict(twinFace, r)) { 
						// The face on the opposite side of trialEdge is also
						// a conflicting (visible) face. Thus, trialEdge does
						// not belong to the outerLoop.
						
						// select the next edge from the current face.
						trialEdge = Polyhedron<T>::getNextEdge(trialEdge); 
						
						//outerLoop.push_back(trialEdge);
						if (debug && ch_debug(r)) {
							std::cout << "trial edge = " << trialEdge << "." << std::endl;
						}
					} else {
	 					// The face on the opposite side of trialEdge is a hidden
 						// face, and trialEdge belongs to the outerLoop.
 						// outerLoop.push_back(trialEdge); // Eureka!
 						initialEdge = trialEdge;
 						if (debug && ch_debug(r)) {
 							std::cout << "Found initialEdge =  " << initialEdge << " to head of outerLoop" << std::endl;
 						}
 					}
					
					// trialEdge = Polyhedron<T>::getNextEdge(trialEdge);
				} while (trialEdge != initialFaceEdge && (! Polyhedron<T>::isDefined(initialEdge)));	
				// Continue iterating until either an edge belonging to the outerLoop
				// has been found, or every edge on every conflicting face has been
				// examined. The latter case should not occur.
			}
			
			//	Ensure that a single edge on the outerLoop has been found. The following
			//	assertion should always be true.
			assert(Polyhedron<T>::isDefined(initialEdge));
		
			// The remainder of the outerLoop is now constructed by navigating
			// counter-clockwise from the edge that was just found.
		
			if (debug && ch_debug(r)) {
					std::cout << "Finding rest of outerLoop." << std::endl;
			}
						
			index_t trialEdge = initialEdge;
			index_t trialFace = Polyhedron<T>::getCurrentFace(trialEdge);
						
			do {
				// Add edge to add to outerloop	
				if (debug && ch_debug(r)) {
					std::cout 	<< "Adding trialEdge = " << trialEdge 
								<< " to the outer loop. " <<  std::endl;
				}
				outerLoop.push_back(trialEdge);	
						
				// Advance search to next vertex.
				trialEdge = Polyhedron<T>::getTwinEdge(trialEdge); 
				index_t edge = trialEdge;
				do {
					// Cycle counter-clockwise around the next vertex of the current trialEdge,
					// until edge falls into a conflicting face.
					edge = Polyhedron<T>::getTwinEdge(Polyhedron<T>::getPrevEdge(edge));
					trialFace = Polyhedron<T>::getCurrentFace(edge);
					if (debug && ch_debug(r)) {
						std::cout << "edge = " << edge 
								  << " belongs to face " << trialFace << std::endl;
					}
					
				    if (! Polyhedron<T>::isDefined(edge)) {
						std::cerr << "Undefined edge implies topological inconsistent polyhedron at r = " 
								  << r << std::endl;
						Polyhedron<T>::describe(std::cerr);
						abort();
					} else if (edge == trialEdge) { // rotated one full revolution without finding
													// the continuation of the horizon.
						std::cerr << "Failed to find outer loop at r = "  << r << std::endl;
						Polyhedron<T>::describe(std::cerr);
						abort();
					}		
				} while (! isaFacePointConflict(trialFace, r));
				
				trialEdge = edge;
			
			} while (trialEdge != initialEdge);
			
			// outerLoop should be complete
					
			if (debug && ch_debug(r)) {
				 std::cout << "done. outerLoop.size() = " << outerLoop.size() << std::endl;
			}
				
			// mark all conflicting faces for deletion. (Note, that the faces are not
			// actually deleted until the polygon is purged.)		
		
			if (debug && ch_debug(r))  std::cout << "Marking faces ";
			for(pos = conFaces.begin(); pos != conFaces.end(); ++pos) {
				if (debug && ch_debug(r)) {
					std::cout << *pos << " "; // the face indices to be marked for deletion.
				}	

				bool status = Polyhedron<T>::deleteFace(*pos);  // Mark the face *pos for deletion.

                if (! status) {
                  std::cerr << "ConvexHull<T>::ConvexHull: ERROR: "
                            << " deleteFace(" << *pos << ") failed at paint " << r << "."
                            << std::endl;
                  Polyhedron<T>::describe(std::cerr);
                  abort();
                }
			}
			if (debug && ch_debug(r)) std::cout << " for deletion." << std::endl;
			
			
			// Ensure that points[r] is in the vertex table.
			index_t vertexIndex;         // the index of the vertex defined by point[r]
			
			if (! Polyhedron<T>::findVertex(points[r], &vertexIndex)) {
			  // point[r] is not in the vertex table, so add it.
			  vertexIndex = Polyhedron<T>::addVertex(points[r]);
			}
		
			// It's now time to replace the deleted faces (those circumscribed by the
			// outerLoop) with triangles subtended by each edge with respect to 
			// point[r].		
			typename std::list<index_t>::iterator loopPos;
			for(loopPos = outerLoop.begin(); loopPos != outerLoop.end(); ++loopPos) {
				// Add the triangle delimited by point[r] and *loopPos to the
				// model.
								
				// Identify the edge and face adjacent to *loopPos:
				index_t adjEdge = Polyhedron<T>::getTwinEdge(*loopPos);
				index_t adjFace = Polyhedron<T>::getCurrentFace(adjEdge);
				index_t loopVertex = Polyhedron<T>::getCurrentVertex(*loopPos);
				
				if (debug && ch_debug(r)) {
					std::cout << "*loopPos   = " << *loopPos << std::endl;
					std::cout << "adjEdge    = " << adjEdge << std::endl;
					std::cout << "adjFace    = " << adjFace << std::endl;
					std::cout << "loopVertex = " << loopVertex << std::endl;
				}
				
				// edgeVector is the displacement to point[r] relative to the
				// current vertex in the loop. It is requrired for the following
				// orthogonality test.
				typename Polyhedron<T>::FaceVertex edgeVector = points[r] 
								- Polyhedron<T>::d_vertexTable[loopVertex].d_v;
								
				double ip = Polyhedron<T>::getNormalVector(adjFace) * edgeVector.reduce();
					
				if (debug && ch_debug(r)) {
					std::cout << "inner product = " << ip << std::endl;
				}
								
				if (std::abs(ip) <= epsilonPerp) {
					// faces are coplanar, so merge them by displacing adjEdge to
					// point[r], and adding a new HalfEdge to close the adjacent
					// face.
					coplanarAdditions++;
					if (debug && ch_debug(r)) {
						std::cout << "Coplanar!" << std::endl;
					}
					
					status = Polyhedron<T>::extendFace(vertexIndex, adjEdge);
											
					if (! status) {
						std::cerr << "ConvexHull<T>: Could not add a new coplanar face corresponding "
							<< "to point with index r = " << r << std::endl;
						Polyhedron<T>::describe(std::cerr);
						abort();
					}
						
					// No need to modify the conflict grqph for adjFace.							
				} else {
					// Faces are not coplanar, so create a new facet.
					// The call to newFace also will add two new edges in the
					// edge table, one from vertexIndex (point[r]) to the
					// endpoint of adjEdge (or origin of *loopPos), and the
					// second from getCurrentVertex(adjEdge) to vertexIndex.
					// addFace also searches for any possible twin edges that
					// might exist, and assigns the twin_ field of the edges
					// accordingly.
					// The value returned is the index of the new face in the
					// face table.
					
					noncoplanarAdditions++;				
					if (debug && ch_debug(r)) { 
						std::cout << "Not coplanar!" << std::endl;
					}
					
					index_t newFace;  // index of the face returned by addFace.
					status = ConvexHull<T>::addFace(	vertexIndex, 
											Polyhedron<T>::getNextVertex(adjEdge), 
											Polyhedron<T>::getCurrentVertex(adjEdge),
											&newFace);
					if (! status) {
						std::cerr << "ConvexHull<T>:: Error:: could not add a new noncoplanar face "
								  << " for point index r = " << r << std::endl;
						Polyhedron<T>::describe(std::cerr);
						abort();
					}
					
					// At this point *loopPos is the twin of adjEdge. Since the former
					// will be deleted, reset the twin of adjEdge to be the edge in the
					// new face that conicides with *loopPos.
					
					// The call to addFace above created three edges. The data member
					// d_outer contains the new edge with the smallest index. The edge
					// coincident with *loopPos is the index that follows d_outer.
					Polyhedron<T>::d_edgeTable[adjEdge].d_twinEdge = 
						Polyhedron<T>::d_faceTable[newFace].d_outer + 1;
					
					// The elements of pointConflicts[newFace] will be a subset
					// of the union of pointConflicts[getCurrentFace(*loopPos)] with
					// pointConflicts[adjFace].
					index_t oldFace = Polyhedron<T>::getCurrentFace(*loopPos);
					
					// std::set<index_t> newPointConflicts;
					
					typename Polyhedron<T>::FaceNormal faceNormal = Polyhedron<T>::getNormalVector(newFace);
					typename Polyhedron<T>::FaceVertex faceVertex = Polyhedron<T>::getFaceVertex(newFace);
					
					if (debug && ch_debug(r)) {
						std::cout << "oldFace = " << oldFace << std::endl;
						std::cout << "adjFace = " << adjFace << std::endl;
						std::cout << "newFace = " << newFace << std::endl;	
						std::cout << "newFace normal = " << faceNormal << std::endl;
						std::cout << "newFace vertex = " << faceVertex << std::endl;
					}
					
					// Inspect the face-point conflict graph
                    if (debug && ch_debug(r)) {
                      this->displayConflictGraph(std::cout, 10);
                    }
					
					typename std::set<index_t>::iterator setPos;
					std::set<index_t> oldFaceConflicts = d_pointConflicts[oldFace];//.getLeftNodeLinkSet(oldFace);
					
					index_t newConflicts = 0;
					for(setPos  = oldFaceConflicts.begin();
						setPos != oldFaceConflicts.end(); ++setPos) {
                      
						if(faceNormal*(points[*setPos] - faceVertex) > epsilonPerp) {
							// Face newFace and point[*setPos] conflict 
							// (both cannot belong to the CH).
							insertFacePointConflict(newFace, *setPos);
							newConflicts++;
						}
					}
					
					std::set<index_t> adjFaceConflicts = d_pointConflicts[adjFace];//.getLeftNodeLinkSet(adjFace);
					for(setPos  = adjFaceConflicts.begin();
						setPos != adjFaceConflicts.end(); ++setPos)
					{
						if( faceNormal*(points[*setPos] - faceVertex) > epsilonPerp) {
							// Face newFace and point[*setPos] conflict 
							// (both cannot belong to the CH).
							insertFacePointConflict(newFace, *setPos);
							newConflicts++;
						}
					}
					
					if (debug && ch_debug(r)) {
						std::cout << "Detected " << newConflicts << " new face-point conflicts." << std::endl;
					}
				} // finished adding a single (coplanar or non-coplanar) face.
				
				// Separate each edge in the loop from the retained edges
				Polyhedron<T>::d_edgeTable[*loopPos].d_twinEdge = Polyhedron<T>::UNDEFINED;
				if (debug && ch_debug(r)) {
					std::cout << std::endl;
				 	this->describe(std::cout);
				}
			}
						
			// Since every face that conflicts with point[r] has been removed from the model,
			// we should delete each vertex in the face-point conflict graph that corresponds
			// to a deleted face, and all links between those nodes and the set of point
			// nodes.

            
            if (debug) {
              if (! isConflictGraphConsistent(points.size()) ) {
                std::cerr << "ConvexHull<T>::ConvexHull: "
                          << "conflict graph is not consistent before purge at r = " << r
                          << std::endl;
                abort();
              }
            }

            if (debug && ch_debug(r)) {
              this->displayConflictGraph(std::cout, 20);
            }

            d_faceConflicts[r].clear(); // r has been added to the CH: it lacks face conflicts.

            purgeConflictGraph();

            // Delete faces in decreasing index order to avoid conflicts that would result
            // from reassigning an index to the last element in the face table.
            
            if (debug && ch_debug(r)) {
              std::cout << "Purging the tables... ";
            }

            Polyhedron<T>::purgeTables();

            if (debug && ch_debug(r)) {
              std::cout << " done." << std::endl;
              this->describe(std::cout);
              this->displayConflictGraph(std::cout, 20);
            }

            if (debug) {
              if (! isConflictGraphConsistent(points.size()) ) {
                std::cerr << "ConvexHull<T>::ConvexHull: "
                          << "conflict graph is not consistent after purge at r = " << r
                          << std::endl;
                abort();
              }


              if (! Polyhedron<T>::checkTables(true)) {
                std::cerr << "ConvexHull::ConvexHull: ERROR: "
                          << "internal tables contain inactive elements after calling purgTables."
                          << std::endl
                          << "Current point is r = " << r
                          << std::endl;
                Polyhedron<T>::describe(std::cerr);
                abort();
              }
            }

			// facePointConflicts.cutNodesAdjacentToRightNode(r);		
            conFaces.clear();
			outerLoop.clear();	// clear the outerLoop in preparation for the next iteration.
            
            if (debug) {
              if (isConsistent) {
                // Test the topological cconstistency of the polyhedron.
                isConsistent = Polyhedron<T>::eulerTest(false);
                if (! isConsistent) {
                  // If the test fails, note the results and cease testing for consistency.
                  // We will handle one bug at time.
                  std::cerr << "Inconsistency detected at r = " << r << std::endl;
                  Polyhedron<T>::eulerTest(true);
                  Polyhedron<T>::describe(std::cerr);
                  abort();
                }
              }
            }
		}	
	}
    std::cout << std::endl; // Terminate the poor man's progress bar

	Polyhedron<T>::eulerTest(true);  // verbose flag is set to true.
}

#endif // __CONVEX_HULL_H__
