// 
//  main.cc
//  polygonTest
//  
//  Created by Robert Snapp on 2007-07-10.
//  Copyright 2007 Robert R. Snapp. All rights reserved.
// 


#include "polygon.h"
#include "ntuple.h"
#include <iostream>
#include <list>

using namespace std;

int main (int argc, char * const argv[]) {
    // insert code here...
	Polygon<double> square;
    vector<Polygon<double>::PolygonVertex> square_pts;
	
	ntuple<double,2> t0(0, 0);
	ntuple<double,2> t1(1, 0);
	ntuple<double,2> t2(1, 1);
	ntuple<double,2> t3(1, 0);
	
	
	square_pts.push_back(t0);
	square_pts.push_back(t1);
	square_pts.push_back(t2);
    square_pts.push_back(t3);
		
    std::cout << "Square:\n";
	square.describe(std::cout);
	
    return 0;
}

