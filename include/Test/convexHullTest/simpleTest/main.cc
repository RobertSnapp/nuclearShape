// 
//  main.cc
//  convexHullTest
//  
//  Created by Robert Snapp on 2007-03-25.
//  Copyright 2007 Robert R. Snapp. All rights reserved.
// 

#include <cstdlib>
#include "polyhedron.h"
#include "ntuple.h"
#include <iostream>
#include <vector>
#include "convexhull.h"

using namespace std;

int const n2 = 10000;  // number of points generated in test #2.
int const n3 = 10000;  // number of points generated in test #3
int const r3[3] = {100, 20, 20}; // Bounds for the rectangular region in test 3.


// Creates a line separator for the output.
inline ostream& divider(ostream& os) {
  for(int i = 0; i < 80; i++) {
	os << "=";
  }
  return os << endl << endl;
}

int main (int argc, char * const argv[]) {

  cout << "First test of convexHull uses the eight vertices of a cube as input:" << endl;

  vector<ntuple<float, 3> > fpoints;
  for(int i = 0; i <= 1; i++) {
	float x = static_cast<float>(i);
	for(int j = 0; j <= 1; j++) {
	  float y = static_cast<float>(j);
	  for(int k = 0; k <= 1; k++) {
		float z = static_cast<float>(k);
		ntuple<float,3> p = ntuple<float,3>(x,y,z);
		fpoints.push_back(p);
	  }
	}
  }

  ConvexHull<float> cube(fpoints);
  cube.describe(cout);
  if (! cube.isConvex()) {
	cout << "ERROR: cube is not convex." << endl;
  }

  divider(cout);

  cout << "Second test computes the convex hull of" << n2 
	   << " random points of type double that lie within a unit cube." << endl;

  cout << "(Here we use rand() with RAND_MAX = " << RAND_MAX << ".)" << endl;
	
  vector<ntuple<double, 3> > dpoints;
  for (int i = 0; i < n2; i++) {
	ntuple<double, 3> x;
	for(int j = 0; j < 3; j++) {
	  x[j] = static_cast<double>(rand())/RAND_MAX;
	}
	dpoints.push_back(x);
  }
	
  ConvexHull<double> hull2(dpoints);
  hull2.describe(cout);
	
  if (hull2.isConvex()) {
	cout << "hull is convex." << endl;
  } else {
	cout << "hull is nonconvex (could be roundoff error)." << endl;
  }

  divider(cout);

  cout << "The third test computes the convex hull of " << n3
	   << " random points of type int that fall within a rectangular volume "
	   << " of dimensions" << endl
	   << "  xmax = " << r3[0] << endl 
	   << "  ymax = " << r3[1] << endl
	   << "  zmax = " << r3[2] << endl;


  vector<ntuple<int, 3> > ipoints;
  for(int i = 0; i < n3; i++) {
	ntuple<int, 3> x;
	for(int j = 0; j < 3; j++) {
	  x[j] = rand() % r3[j];
	}
	ipoints.push_back(x);
  }

  ConvexHull<int> hull3(ipoints);
  hull3.describe(cout);

  return 0;
}
