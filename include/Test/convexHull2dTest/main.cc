// 
//  main.cc
//  convexHull2dTest
//  
//  Created by Robert Snapp on 2007-07-11.
//  Copyright 2007 Robert R. Snapp. All rights reserved.
// 

#include <cstdlib>
#include "polygon.h"
#include "ntuple.h"
#include <iostream>
#include <vector>
#include "convexHull2d.h"

using namespace std;

int main (int argc, char * const argv[]) {
	
	vector<ntuple<float, 2> > points;
	
	ntuple<float,2> c0(0,0);
	ntuple<float,2> c1(1,0);
	ntuple<float,2> c2(2,1);
	ntuple<float,2> c3(0,2);
	ntuple<float,2> c4(1,1);
	ntuple<float,2> c5(0,3);
	ntuple<float,2> c6(2,2);
	ntuple<float,2> c7(3,0);
	
	points.push_back(c0);
	points.push_back(c1);
	points.push_back(c2);
	points.push_back(c3);
	points.push_back(c4);
	points.push_back(c5);
	points.push_back(c6);
	points.push_back(c7);
	
	ConvexHull2d<float> factory(points);
	
	if (factory.isConvex(true)) {
		cout << "first polygon is convex." << endl;
	} else {
		cout << "first polygon is nonconvex." << endl;
	}
	

	factory.describe(cout);
	
	cout << "RAND_MAX = " << RAND_MAX << endl;
	
	vector<ntuple<double, 2> > dpoints;
	for (int i = 0; i < 10000; i++) {
		ntuple<double, 2> x;
		for(int j = 0; j < 2; j++) {
			x[j] = static_cast<double>(rand())/RAND_MAX;
		}
		dpoints.push_back(x);
	}
	
	ConvexHull2d<double> f2(dpoints);
	f2.describe(cout);
	
	if (f2.isConvex(true)) {
		cout << "Polygon f2 is convex." << endl;
	} else {
		cout << "Polygon f2 is nonconvex." << endl;
	}
	
    return 0;
}
