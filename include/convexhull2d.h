#ifndef __CONVEX_HULL_2D_H__
#define __CONVEX_HULL_2D_H__

//
// convexHull2d.h
//
// Created by Robert R. Snapp on 2007-07-10.
// Copyright (C) Robert R. Snapp. All rights reserved.
//

#include <vector>
#include <set>
#include <iostream>
#include <ostream>
#include "polygon.h"

// ch_debug(r) facilitates testing and debuging conditioned on the index of the current point under
// consideration. If it returns true, then a lot of information about the current geometric model
// will be displayed in the console.
inline bool ch_debug(int r) {
#pragma unused(r)
  return false;
}

template <typename T>
class ConvexHull2d : public Polygon<T> {
public: 
  ConvexHull2d();
  ConvexHull2d(std::vector<typename Polygon<T>::PolygonVertex> &points);
  void compute(std::vector<typename Polygon<T>::PolygonVertex> &points);
 
protected:
  bool presortConvexVector( std::vector<typename Polygon<T>::PolygonVertex> &points, 
                            double conditionFactor);
  void insertEdgePointConflict(index_t e, index_t p);
  bool isConflictGraphConsistent(index_t npts) const;
  bool isanEdgePointConflict(index_t f, index_t p) const;
  bool addEdge(index_t v1, index_t prev, index_t next, index_t* newEdge);
  bool addEdge(Point2d &v, index_t prev, index_t next, index_t* newEdge);
  bool addEdge(typename Polygon<T>::PolygonVertex &v, index_t prev, index_t next, index_t* newEdge);
  bool purgeConflictGraph(void);
  void displayConflictGraph(std::ostream &os, index_t maxIndex) const;

  // The i-th element of edgeConflicts should return the list of edge-indicies
  // that conflict with point[i].
  std::vector<std::set<index_t> > edgeConflicts_;
  
  // The j-th element of pointConflicts should return the list of point-indices
  // that conflict with edge[j]
  std::vector<std::set<index_t> > pointConflicts_;
};


template<typename T>
void ConvexHull2d<T>::insertEdgePointConflict(index_t e, index_t p) {
  edgeConflicts_[p].insert(e);
  pointConflicts_[e].insert(p);
}

template<typename T>
bool ConvexHull2d<T>::isConflictGraphConsistent(index_t npts) const {
  if (Polygon<T>::edges() != pointConflicts_.size()) {
    std::cerr << "ConvexHull2d<T>::isConflictGraphConsistent: "
              << "edges() = " << Polygon<T>::edges() 
              << "but pointConflicts_.size() = " << pointConflicts_.size()
              << std::endl;
    return false;
  }

  if (npts != edgeConflicts_.size()) {
    std::cerr << "ConvexHul2dl<T>::isConflictGraphConsistent: "
              << "npts = " << npts 
              << "but edgeConflicts_.size() = " << edgeConflicts_.size()
              << std::endl;
    return false;
  }

  typename std::set<index_t>::iterator p;
  for (index_t i = 0; i < Polygon<T>::edges(); i++) {
    for(p = pointConflicts_[i].begin(); p != pointConflicts_[i].end(); ++p) {
      if (edgeConflicts_[*p].find(i) == edgeConflicts_[*p].end()) {
        std::cerr << "ConvexHull2d<T>::isConflictGraphConsistent: "
                  << "edge inconsistency between edge i = " << i
                  << " and point *p = " << *p
                  << std::endl;
        return false;
      }
    }
  }

  typename std::set<index_t>::iterator e;
  for (index_t i = 0; i < npts; i++) {
    for(e = edgeConflicts_[i].begin(); e != edgeConflicts_[i].end(); ++e) {
      if (pointConflicts_[*e].find(i) == pointConflicts_[*e].end()) {
        std::cerr << "ConvexHull2d<T>::isConflictGraphConsistent: "
                  << "edge inconsistency between point i = " << i
                  << " and edge *e = " << *e
                  << std::endl;
        return false;
      }
    }
  }

  return true;
}
  
template<typename T>
bool ConvexHull2d<T>::isanEdgePointConflict(index_t f, index_t p) const {
  return (edgeConflicts_[p].find(f) != edgeConflicts_[p].end());
}

template<typename T>
bool ConvexHull2d<T>::addEdge(index_t v1, index_t prev, index_t next, index_t* edgeIndex) {
  bool status;
  pointConflicts_.push_back(std::set<int>());

  status = Polygon<T>::addEdge(v1, prev, next, edgeIndex);
  if (!status) {
    std::cerr << "ConvexHull2d<T>::addEdge(" << v1 << ", " << prev << ", " << next << ",...): "
              << "call to Polygon<T>::addEdge(...) failed."
              << std::endl;
    return false;
  }

  if (Polygon<T>::edges() != pointConflicts_.size()) {
    std::cerr << "ConvexHull2d<T>::addEdge(" << v1 << ", " << prev << ", " << next << ",...): "
              << "pointConflicts_.size() = " << pointConflicts_.size()
              << " but edges() = " << Polygon<T>::edges() 
              << std::endl;
    return false;
  }

  return status;
}

template<typename T>
bool ConvexHull2d<T>::addEdge(typename Polygon<T>::PolygonVertex &v1, index_t prev, index_t next, index_t* edgeIndex) {
  bool status;
  pointConflicts_.push_back(std::set<index_t>());

  status = Polygon<T>::addEdge(v1, prev, next, edgeIndex);
  if (!status) {
    std::cerr << "ConvexHull2d<T>::addEdge(" << v1 << ", " << prev << ", " << next << ",...): "
              << "call to Polygon<T>::addEdge(...) failed."
              << std::endl;
    return false;
  }

  if (Polygon<T>::edges() != pointConflicts_.size()) {
    std::cerr << "ConvexHull2d<T>::addEdge(" << v1 << ", " <<  prev << ", " << next << ",...): "
              << "pointConflicts_.size() = " << pointConflicts_.size()
              << " but edges() = " << Polygon<T>::edges() 
              << std::endl;
    return false;
  }

  return status;
}

template<typename T>
bool ConvexHull2d<T>::addEdge(Point2d &v1, index_t prev, index_t next, index_t* edgeIndex) {
  return ConvexHull2d<T>::addEdge(ntuple<T,2>(v1), prev, next, edgeIndex);
}

// Function purgeConflictGraph modifies the edge-point conflict graph to account for
// the deletions of the edges that are members of the set deadEdges_.  The boolean
// value true is returned if the opeartion was performed successfully.
template<typename T>
bool ConvexHull2d<T>::purgeConflictGraph(void) {
  bool status = true;

  typename Polygon<T>::IntSet::iterator e;
  for (e = Polygon<T>::deadEdges_.begin(); e != Polygon<T>::deadEdges_.end(); ++e) {
    // e indicated an edge that is being deleted.
    index_t eLast = pointConflicts_.size() - 1;

    // Remove all references to *e from the edgeConflict_ link sets: indexed by members of 
    // pointConflicts_[*e].
    typename std::set<index_t>::iterator p;
    for(p = pointConflicts_[*e].begin(); p != pointConflicts_[*e].end(); ++p) {
      edgeConflicts_[*p].erase(*e);
    }

    if (*e != eLast) {
      // copy the last item in the pointConflicts table, overwriting the record pointConflicts_[*e].
      pointConflicts_[*e] = pointConflicts_[eLast];
     
      // Change all references of edge eLast, to references of edge *e.
      for(p = pointConflicts_[eLast].begin(); p != pointConflicts_[eLast].end(); ++p) {
        if (edgeConflicts_[*p].erase(eLast) > 0) {
          edgeConflicts_[*p].insert(*e);
        } else {
          status = false;
        }
      }
    }
    pointConflicts_.pop_back();
  }

  return status;
}

template<typename T>
void ConvexHull2d<T>::displayConflictGraph(std::ostream &os, index_t maxIndex) const {
  os << "Conflict Graph: (" 
     << pointConflicts_.size() << " edge nodes, "
     << edgeConflicts_.size()  << " point nodes)"
     << std::endl;
  
  for(index_t i = 0; i < std::min(pointConflicts_.size(), maxIndex); i++) {
    os << "L" << i << ": ";
    index_t j = 0;
    std::set<index_t>::iterator p;
    for(p = pointConflicts_[i].begin(); p != pointConflicts_[i].end() && j < maxIndex + 1; ++p, ++j) {
      os << *p << ", ";
      if (j == maxIndex) {
        os << "...";
      } 
    }
    os << std::endl;
  }

  os << "-----------------"
     << std::endl;

  for(index_t i = 0; i < std::min(edgeConflicts_.size(), maxIndex); i++) {
    os << "R" << i << ": ";
    index_t j = 0;
    std::set<index_t>::iterator e;
    for(e = edgeConflicts_[i].begin(); e != edgeConflicts_[i].end() && j < maxIndex + 1; ++e, ++j) {
      os << *e << ", ";
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
bool ConvexHull2d<T>::presortConvexVector(
                                          std::vector<typename Polygon<T>::PolygonVertex> &points, 
                                          double conditionFactor) {
  assert(0 <= conditionFactor && conditionFactor < 1.0);
  index_t size = points.size();
	
  // Randomize the vector points:
  for(index_t count = 0; count < size; count++) {
    index_t p = (std::rand() % size);
    index_t q = (std::rand() % size);
    std::swap(points[p], points[q]);
  }
	
  for(index_t i0 = 0; i0 < size - 2; ++i0)
    for(index_t i1 = i0 + 1; i1 < size - 1; ++i1)
      for(index_t i2 = i1 + 1; i2 < size; ++i2) {
        typename Polygon<T>::PolygonVertex v1 = points[i1] - points[i0]; // sides of the 
        typename Polygon<T>::PolygonVertex v2 = points[i2] - points[i0]; // parallelogram
        					
        double bound = conditionFactor*v1.norm()*v2.norm();
        if (Polygon<T>::debug && ch_debug(0)) {
          std::cerr << "bound = " << bound;
        }
					
        double det   = v1 * perp(v2); // oriented area
        if (Polygon<T>::debug && ch_debug(0)) {
          std::cerr << ", parallelogram area = " << det << std::endl;
        }
					
        if (std::abs(det) > bound) {
          if (Polygon<T>::debug && ch_debug(0)) {
            std::cerr << "i0 = " << i0
                      << ", i1 = " << i1 
                      << ", i2 = " << i2 
                      << std::endl;
          }
          std::swap(points[0], points[i0]);
          std::swap(points[1], points[i1]);
          std::swap(points[2], points[i2]);
						
          if (det < 0) {
            std::swap(points[1], points[2]);
          }
						
          v1 = points[1] - points[0]; // sides of the parallelogram
          v2 = points[2] - points[0];
						
          assert(v1 * perp(v2) > 0);
          return true; // Exits with success
        }
      }
  return false; // Exits with failure.
}

template <typename T>
ConvexHull2d<T>::ConvexHull2d() : Polygon<T>::Polygon() {}

template <typename T>
ConvexHull2d<T>::ConvexHull2d(std::vector<typename Polygon<T>::PolygonVertex> &points) : Polygon<T>::Polygon() {
  compute(points);
}

template <typename T>
void ConvexHull2d<T>::compute(std::vector<typename Polygon<T>::PolygonVertex> &points) {
  this->clear();

  bool isConsistent = true;  // This flag is used for debugging. It should be set to false
  // if a topological inconsistency is detected by eulerTest.
  double conditionFactor = 0.20;

  if (points.size() < 3) {
    std::cerr << "ConvexHull2d<T>::ConvexHull2d(std::vector<PolygonVertex>): "
              << "Number of input points must exceed 2." << std::endl;
    exit(1);
  }
	
  // presortConvexVector reorders the vector points so that the first three elements form a
  // non-degenerate triangle.
  bool status = presortConvexVector(points, conditionFactor);
  assert(status);
  
  // Polygon<T>::clear();
  // Construct the initial triangle and verify that edge orientations are correct.
  Polygon<T>::createTriangle(points[0], points[1], points[2]);

  assert(Polygon<T>::edges()    == 3);
  assert(Polygon<T>::vertices() == 3);
	
  if (Polygon<T>::debug && ch_debug(0)) {
    this->describe(std::cout);
    std::cout << std::endl;
  }
	
  // The conflict graph is a bipartite graph that consists of two types of
  // nodes. On set of nodes, pointConflicts, contains one node per polygon
  // edge. Thus, pointConflicts[i] refers to the subset of input points that 
  // conflict with edge i. The other set of nodes, edgeConflicts_, contains
  // one node per ntuple in points. Thus, edgeConflicts_[j] refers to the 
  // set of edges that conflict with points[j].
	
  // Initialize the conflict graph:
  if (Polygon<T>::debug && ch_debug(0)) {
    std::cout << "Initializing the conflict graph ...";
  }

  edgeConflicts_.resize(points.size(), std::set<index_t>());
  pointConflicts_.resize(3, std::set<index_t>());
   
  int conflicts = 0;
  for(index_t e = 0; e < 3; ++e) {
    for(index_t p = 3; p < points.size(); ++p) {
      if(Polygon<T>::isVisible(e, points[p])) {
        // Edge e and point p conflict (both cannot belong to the CH).
        conflicts++;
        ConvexHull2d<T>::insertEdgePointConflict(e,p);
      }
    }
  }
		
  if (Polygon<T>::debug && ch_debug(0)) {
    std::cout << "done: " <<  conflicts << " edges.\n";
    if (isConflictGraphConsistent(points.size())) {
      std::cout << "edge-point conflict graph passes initial symmetry test." << std::endl;
    } else {
      std::cout << "edge-point conflict graph fails initial symmetry test." << std::endl;
    }
  }
	
  // Incrementally, construct the convex hull (CH).
  for(index_t r = 3; r < points.size(); ++r) {
    // remove all edges that are in conflict with points[r]
	
    if (Polygon<T>::debug) {
      // Logarithmic "poor man's" progress bar.
      if (r %    10 == 0 && r < 250) std::cout << "." << std::flush;
      if (r %   100 == 0 && r < 1000) std::cout << ","<< std::flush; 
      if (r %  1000 == 0 && r < 10000 && r >= 1000) std::cout << ":"<< std::flush;
      if (r % 10000 == 0 && r >= 10000) std::cout << "!"<< std::flush;
    }

    if (edgeConflicts_[r].size()  > 0) {
      // point[r] is outside the current convex polygon and must
      // be added to the convex hull.
      if (Polygon<T>::debug && ch_debug(r)) {
        std::cout  << std::endl
                   << "Adding points[" << r << "]:" << std::endl
                   << "points[" << r << "] = " << points[r] 
                   << std::endl
                   << std::endl;
        Polygon<T>::describe(std::cout);
        std::cout << std::endl;
      }
		
      // Identify the first and last edges that delemits the set of conflicting
      // edges, i.e., those that are visible from point[r].

      if (Polygon<T>::debug && ch_debug(r)) {
        std::cout << "Constructing the outer loop for points[" << r << "] ...\n";
      }
			
      // The set of edges that conflict with points[r] should be contiguous. Let
      // headEdge denote the first (leftmost) conflicting edge, and let lastEdge
      // denote the last (rightmost) conflicting edge. We will also let tailEdge
      // denote the successor to lastEdge. The following statements are designed
      // to determine the values of headEdge, lastEdge, and tailEdge:
			
      std::set<index_t>::iterator pos = edgeConflicts_[r].begin();
            
      index_t headEdge = *pos;  // A hypothetical initial value.
      index_t prevEdge;
      do {
        prevEdge = Polygon<T>::getPrevEdge(headEdge);
        if (edgeConflicts_[r].find(prevEdge) != edgeConflicts_[r].end()) {
          // prevEdge is also in edgeConflicts.
          headEdge = prevEdge;
        }

        // If prevEdge does not conflict with points[r], then we stop
        // itorating.
      } while (headEdge == prevEdge);
      // At this point headEdge is determined. Now to find lastEdge:

      if (Polygon<T>::debug && ch_debug(r)) {
        std::cerr << "Head edge = " << headEdge << std::endl;
      }
            
      index_t lastEdge = *pos; // A hypothetical initial value
      index_t tailEdge;
      do {
        tailEdge = Polygon<T>::getNextEdge(lastEdge);
        if (edgeConflicts_[r].find(tailEdge) != edgeConflicts_[r].end()) {
          lastEdge = tailEdge;
        }
      } while (lastEdge == tailEdge);
      // At this point, lastEdge and tailEdge are determined.

      if (Polygon<T>::debug && ch_debug(r)) {
        std::cerr << "Last edge = " << lastEdge << std::endl
                  << "Tail edge = " << tailEdge << std::endl;
      }

      // To update the two-dimensional convex hull, delete edges from firstConflict->next through lastConflict.

      // Since edge headEdge will become incident to points[r], this conflict should be removed
      // from the conflict graph.
      edgeConflicts_[r].erase(headEdge);
      pointConflicts_[headEdge].erase(r);

      // mark the remaining conflicting edges for deletion. (Note, that the edges are not
      // actually deleted until the tables are purged.)	
      if (headEdge != lastEdge) {

        // Update the initialEdge_ pointer to ensure that it references a surviving edge.
        if (edgeConflicts_[r].find(Polygon<T>::initialEdge_) != edgeConflicts_[r].end()) {
          Polygon<T>::initialEdge_ = headEdge;
        }

        if (Polygon<T>::debug && ch_debug(r))  std::cout << "Marking edges ";
        for (pos = edgeConflicts_[r].begin(); pos != edgeConflicts_[r].end(); ++pos) {
          if (Polygon<T>::debug && ch_debug(r)) {
            std::cout << *pos << " ";
          }

          status = Polygon<T>::deleteEdge(*pos);
          if (! status) {
            std::cerr << "ConvexHull2d<T>::ConvexHull2d: ERROR: "
                      << " deleteEdge(" << *pos << ") failed at paint " << r << "."
                      << std::endl;
            Polygon<T>::describe(std::cerr);
            abort();
          }
        }
        if (Polygon<T>::debug && ch_debug(r)) std::cout << " for deletion." << std::endl;
      }	
	
      // Now we create a new edge that will join points[r] with the tailEdge. Also, add a new
      // edge vertex in the edge-point conflict graph.
      index_t newEdge;
      status = ConvexHull2d<T>::addEdge(points[r], headEdge, tailEdge , &newEdge);

      if (! status) {
        std::cerr << "ConvexHull2d<T>::ConvexHull2d(...): "
                  << "Could not add edge for point r = " << r
                  << std::endl;
        abort();
      }

      if (Polygon<T>::debug && ch_debug(r)) {
        std::cerr << "Updated tables:"
                  << std::endl << std::endl;
        this->describe(std::cerr);
        std::cerr << std::endl;
      }

            
      // Inspect the edge-point conflict graph
      if (Polygon<T>::debug && ch_debug(r)) {
        std::cerr << "Previous conflict graph: "
                  << std::endl;
        this->displayConflictGraph(std::cout, 10);
        std::cerr << std::endl;
      }

      // Update the edge-point conflict graph:
      // First, construct conflicts for newEdge
      for (pos = pointConflicts_[headEdge].begin(); pos != pointConflicts_[headEdge].end(); ++pos) {
        if (Polygon<T>::isVisible(newEdge, points[*pos])) {
          ConvexHull2d<T>::insertEdgePointConflict(newEdge, *pos);
        }
      }
           
      // Second modify the point conflicts for headEdge. (A copy of pointConflicts_[headEdge]
      // is created so that the new modifications don't corrupt the set traversal.)
      std::set<index_t> conflictSet = pointConflicts_[headEdge];
      for (pos = conflictSet.begin(); pos != conflictSet.end(); ++pos) {
        if (! Polygon<T>::isVisible(headEdge, points[*pos])) {
          pointConflicts_[headEdge].erase(*pos);
          edgeConflicts_[*pos].erase(headEdge);
        }
      }

      index_t prev = Polygon<T>::getPrevEdge(headEdge);
      for (pos = pointConflicts_[prev].begin(); pos != pointConflicts_[prev].end(); ++pos) {
        if (Polygon<T>::isVisible(headEdge, points[*pos])) {
          ConvexHull2d<T>::insertEdgePointConflict(headEdge, *pos);
        }
      }

      for (pos = pointConflicts_[lastEdge].begin(); pos != pointConflicts_[lastEdge].end(); ++pos) {
        if (Polygon<T>::isVisible(newEdge, points[*pos])) {
          ConvexHull2d<T>::insertEdgePointConflict(newEdge, *pos);
        }
      }

       for (pos = pointConflicts_[tailEdge].begin(); pos != pointConflicts_[tailEdge].end(); ++pos) {
        if (Polygon<T>::isVisible(newEdge, points[*pos])) {
          ConvexHull2d<T>::insertEdgePointConflict(newEdge, *pos);
        }
      }


               						
      // Since every edge that conflicts with point[r] has been removed from the model,
      // we should delete each vertex in the edge-point conflict graph that corresponds
      // to a deleted edge, and all links between those nodes and the set of point
      // nodes.
 
      if (Polygon<T>::debug && ch_debug(r)) {
        this->displayConflictGraph(std::cout, 30);
      }

      if (Polygon<T>::debug) {
        if (! isConflictGraphConsistent(points.size()) ) {
          std::cerr << "ConvexHull2d<T>::ConvexHull2d(...): "
                    << "conflict graph is not consistent before purge at r = " << r
                    << std::endl;
          abort();
        }
      }

  
      edgeConflicts_[r].clear(); // r has been added to the CH: it lacks edge conflicts.

      purgeConflictGraph();

      // Delete edges in decreasing index order to avoid conflicts that would result
      // from reassigning an index to the last element in the edge table.
            
      if (Polygon<T>::debug && ch_debug(r)) {
        std::cout << "Purging the tables... ";
      }

      Polygon<T>::purgeTables();

      if (Polygon<T>::debug && ch_debug(r)) {
        std::cout << " done." << std::endl;
        this->describe(std::cout);
        this->displayConflictGraph(std::cout, 20);
      }

      if (Polygon<T>::debug) {
        if (! isConflictGraphConsistent(points.size()) ) {
          std::cerr << "ConvexHull<T>::ConvexHull: "
                    << "conflict graph is not consistent after purge at r = " << r
                    << std::endl;
          abort();
        }


        if (! Polygon<T>::checkTables()) {
          std::cerr << "ConvexHull::ConvexHull: ERROR: "
                    << "internal tables contain inactive elements after calling purgeTables."
                    << std::endl
                    << "Current point is r = " << r
                    << std::endl;
          Polygon<T>::describe(std::cerr);
          abort();
        }
      }
            
      if (Polygon<T>::debug && ch_debug(r)) {
        if (Polygon<T>::isConvex(true)) {
          std::cerr << "After r = " << r << " polygon is convex\n";
        } else {
           std::cerr << "After r = " << r << " polygon is NOT convex\n";
        }

        if (isConsistent) {
          // Test the topological cconstistency of the polyhedron.
          isConsistent = Polygon<T>::topologyTest(false);
          if (! isConsistent) {
            // If the test fails, note the results and cease testing for consistency.
            // We will handle one bug at time.
            std::cerr << "Inconsistency detected at r = " << r << std::endl;
            Polygon<T>::topologyTest(true);
            Polygon<T>::describe(std::cerr);
            abort();
          }
        }
      }
    }	
  }

  edgeConflicts_.clear();
  pointConflicts_.clear();

  if (Polygon<T>::debug) {
    Polygon<T>::topologyTest(true);  // verbose flag is set to true.
  }
}

#endif // __CONVEX_HULL_2D_H__
