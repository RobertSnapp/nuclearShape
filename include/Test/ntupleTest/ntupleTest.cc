//
//  ntupleTest.cc
//
//  Created by Robert Snapp on 2007-02-16.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include "ntuple.h"
#include <iostream>
using namespace std;

// errorLog will be used to print a message to cerr everytime a scheduled test fails.
// The second argument should indicate the location (i.e., line number) of the test.
inline void errorLog(char* const testName, int line) {
  std::cerr << "ERROR: " << testName 
			<< " test failed above line " 
			<< line 
			<< " in file " 
			<< __FILE__ 
			<< std::endl;
}

// The following macro uses the C++ preprocessor to determine the line number, as the macro __LINE__ is
// evaluted *after* the macro substitution.
#define ERROR_LOG(testname) errorLog(testname, __LINE__)

// Some of the tests involve validating an equality between pairs of floating point numbers.
// Since finite-precision arithmetic complicates this, the following function validates
// the equality of two double precision variables within a preprescribed epsilon.
const double epsilon = 1.0E-08;
inline bool withinEpsilon(double x, double y) {
  return abs(x-y) < epsilon;
}


int main(int argc, char **argv) {
	cerr << "Testing template ntuple:" << endl;
	int failures  = 0; // The total number of failed tests.
	int successes = 0; // The total number of successful tests.
	bool testStatus;   // flag use for multiple-part tests: true indicates success; false, failure.
	
	//////////////////////////////////////////////////////////////////////
    // Test Constructors that initialize from a C++ array.

	int    x6data[] = {1, 2, 3, 4, 5, 6};
	unsigned char ucdata[] = {1, 2, 3, 4, 5, 6};

	float  y6data[6];
	double z6data[6];

	for (int i = 0; i < 6; i++) {
		y6data[i] = 1.0/static_cast<double>(x6data[i]);
		z6data[i] = (i % 2 == 1 ? -1.0 : 1.0)/static_cast<double>(x6data[i]);
	}

	ntuple<int, 6>    ivec(x6data);
	ntuple<float, 6>  fvec(y6data);
	ntuple<double, 6> dvec(z6data);

	// (The complete test is implicit in what follows.)
	
	//////////////////////////////////////////////////////////////////////////////////////
	// Test Constructors that initialize from provided values.

	ntuple<int, 3> zero3(0,0,0);
	ntuple<int, 3> ihat(1,0,0);
	ntuple<int, 3> jhat(0,1,0);
	ntuple<int, 3> khat(0,0,1);

	// (The complete test is implicit in what follows.)

	///////////////////////////////////////////////////////////
	// Test constructors that provide incomplete arguments.
	ntuple<long, 6> fourArguments(0, 1, 2, 3); // The last two components should be zero.
	testStatus = true; // initialize 

	for(int i = 0; i < 6; i++) {
	  testStatus &= (fourArguments[i] == (i < 4? i : 0));
	}

	if (! testStatus) {
	  ERROR_LOG("Incomplete constuctor");
	  failures++;
	} else {
	  successes++;
	}

	//////////////////////////////////////////////////////////
	// Test assignment operator=

	ntuple<int, 6> ivecCopy;
	ivecCopy = ivec; // homogeneous int assignment

	ntuple<double, 6> dvecCopy;
	dvecCopy = dvec;

	testStatus = true; // initialize error flag.

	for (int i = 0; i < 6; i++) {
	  testStatus &= (ivecCopy[i] == ivec[i]);
	}

	if (! testStatus) {
	  ERROR_LOG("Assignment");
	  failures++;
	} else {
	  successes++;
	}

	/////////////////////////////////////////////////////////////
	// Test member function clear()
	ivecCopy.clear();
	testStatus = true;
	for(int i = 0; i < 6; i++) {
	  testStatus &= (ivecCopy[i] == 0);
	}

	if (! testStatus) {
	  ERROR_LOG("clear() ");
	  failures++;
	} else {
	  successes++;
	}
	/////////////////////////////////////////////////////////////
	// Test arithmetic and conditional operations.
	ntuple<int, 3> ijk = ihat + jhat + khat;
	
	cerr << "Testing operator<< :" << endl;
	cout << "ijk = " << ihat << " + " << jhat << " + " << khat << " = " << ijk << endl;
	testStatus = (ijk == ntuple<int, 3>(1,1,1));
	if (! testStatus) {
	  ERROR_LOG("ijk == ntuple<int, 3>(1,1,1)");
	  failures++;
	} else {
	  successes++;
	}
	
	testStatus = ! (ijk != ntuple<int, 3>(1,1,1));
	if (! testStatus) {
	  ERROR_LOG("! (ijk != ntuple<int, 3>(1,1,1))");
	  failures++;
	} else {
	  successes++;
	}
	
	// Test scalar-ntuple multiplication:
	ntuple<int, 3> three_ijk = 3*ijk;	
	cout << "three_ijk = 3*ijk = " << three_ijk << endl;
	testStatus = (three_ijk == ntuple<int, 3>(3,3,3));

	if (! testStatus) {
	  ERROR_LOG("three_ijk == ntuple<int, 3>(3,3,3)");
	  failures++;
	} else {
	  successes++;	
	}
	
	// Test ntuple-scalar multiplication
	ntuple<int, 3> three_ijk_alt = ijk*3;
	cout << "three_ijk_alt = ijk*3 = " << three_ijk_alt << endl;

	testStatus = (three_ijk_alt == ntuple<int, 3>(3,3,3));
	if (! testStatus) {
	  ERROR_LOG("three_ijk_alt == ntuple<int, 3>(3,3,3)");
	  failures++;
	} else {
	  successes++;	
	}
	
	// zero3 is declared above to be (0,0,0).
	testStatus = (three_ijk.clear() == zero3);
	if (! testStatus) {
	  ERROR_LOG("three_ijk.clear() == zero3");
	  failures++;	
	} else {
	  successes++;
	}
	
	// testing operator/
	testStatus = ((2*ijk + ijk)/3 == ijk);
	if (! testStatus) {
	  ERROR_LOG("operator/");
	  failures++;
	} else {
	  successes++;
	}
	
	// testing operator+=
	three_ijk = 3*ijk;
	three_ijk += ijk;
	ntuple<int, 3> four_ijk = three_ijk;
	testStatus = (four_ijk == ntuple<int, 3>(4,4,4));
	if (! testStatus) {
	  ERROR_LOG("operator+=");
	  failures++;
	} else {
	  successes++;
	}
	
	// testing operator-=
	three_ijk = ntuple<int,3>(3,3,3);
	four_ijk -= ijk;
	testStatus = (four_ijk == three_ijk);
	if (! testStatus) {
	  ERROR_LOG("operator-=");
	  failures++;
	} else {
	  successes++;
	}
	
	// testing innerproduct (with integer ntuple)
	int ihat_dot_jhat = ihat.innerProduct(jhat);
	int ihat_dot_khat = ihat.innerProduct(khat);
	int jhat_dot_khat = jhat.innerProduct(khat);
	int ihat_dot_ihat = ihat.innerProduct(ihat);
	int jhat_dot_jhat = jhat.innerProduct(jhat);
	int khat_dot_khat = khat.innerProduct(khat);
	
	if (ihat_dot_jhat == 0 && 
		ihat_dot_khat == 0 &&
		jhat_dot_khat == 0 &&
		ihat_dot_ihat == 1 &&
		jhat_dot_jhat == 1 &&
		khat_dot_khat == 1   ) {
		cout << "Passed integer inner product test." << endl;
		successes++;
	} else {
		cout << "ERROR: Failed integer inner product test." << endl;
		failures++;
	}
	cout << endl;
	
	// testing innerproduct (with float ntuple)
	ntuple<float,3> xhat(1.0, 0.0, 0.0);
	ntuple<float,3> yhat(0.0, 1.0, 0.0);
	ntuple<float,3> zhat(0.0, 0.0, 1.0);
	float epsilon = 1.0E-05;
	
	float xhat_dot_yhat = xhat.innerProduct(yhat);
	float xhat_dot_zhat = xhat.innerProduct(zhat);
	float yhat_dot_zhat = yhat.innerProduct(zhat);
	float xhat_dot_xhat = xhat.innerProduct(xhat);
	float yhat_dot_yhat = yhat.innerProduct(yhat);
	float zhat_dot_zhat = zhat.innerProduct(zhat);
	
	if (abs(xhat_dot_yhat) <= epsilon && 
		abs(xhat_dot_zhat) <= epsilon &&
		abs(yhat_dot_zhat) <= epsilon &&
		abs(xhat_dot_xhat - 1) <= epsilon &&
		abs(yhat_dot_yhat - 1) <= epsilon &&
		abs(zhat_dot_zhat - 1) <= epsilon    ) {
		cout << "Passed float inner product test." << endl;
		successes++;
	} else {
		cout << "ERROR: Failed float inner product test." << endl;
		cout << xhat << "*" << yhat << " = " << xhat_dot_yhat << endl;
		cout << xhat << "*" << zhat << " = " << xhat_dot_zhat << endl;
		cout << yhat << "*" << zhat << " = " << yhat_dot_zhat << endl;
		cout << xhat << "*" << xhat << " = " << xhat_dot_xhat << endl;
		cout << yhat << "*" << yhat << " = " << yhat_dot_yhat << endl;
		cout << zhat << "*" << zhat << " = " << zhat_dot_zhat << endl;
		
		failures++;
	}
	cout << endl;
	
	// testing innerproduct (with int and float ntuple)
	float ihat_dot_yhat = ihat.innerProduct(yhat);
	float xhat_dot_khat = xhat.innerProduct(khat);
	float jhat_dot_zhat = jhat.innerProduct(zhat);
	float xhat_dot_ihat = xhat.innerProduct(ihat);
	float jhat_dot_yhat = jhat.innerProduct(yhat);
	float zhat_dot_khat = zhat.innerProduct(khat);
	
	if (abs(ihat_dot_yhat) <= epsilon && 
		abs(xhat_dot_khat) <= epsilon &&
		abs(jhat_dot_zhat) <= epsilon &&
		abs(xhat_dot_ihat - 1) <= epsilon &&
		abs(jhat_dot_yhat - 1) <= epsilon &&
		abs(zhat_dot_khat - 1) <= epsilon    ) {
		cout << "Passed float-int, int-float inner product test." << endl;
		successes++;
	} else {
		cout << "ERROR: Failed float-int, int-float inner product test." << endl;
		cout << ihat << "*" << yhat << " = " << ihat_dot_yhat << endl;
		cout << xhat << "*" << khat << " = " << xhat_dot_khat << endl;
		cout << jhat << "*" << zhat << " = " << jhat_dot_zhat << endl;
		cout << xhat << "*" << ihat << " = " << xhat_dot_ihat << endl;
		cout << jhat << "*" << yhat << " = " << jhat_dot_yhat << endl;
		cout << zhat << "*" << khat << " = " << zhat_dot_khat << endl;
		
		failures++;
	}
	cout << endl;
	
	// testing innerproduct operator*
	if (static_cast<int>(ihat*jhat) == 0 && 
		static_cast<int>(ihat*khat) == 0 &&
		static_cast<int>(jhat*khat) == 0 &&
		static_cast<int>(ihat*ihat) == 1 &&
		static_cast<int>(jhat*jhat) == 1 &&
		static_cast<int>(khat*khat) == 1 &&
		abs(xhat*jhat) < epsilon         &&
		abs(ihat*yhat) < epsilon         &&
		abs(xhat*yhat) < epsilon           ) {
		cout << "Passed  innerproduct operator* test." << endl;
		successes++;
	} else {
		cout << "ERROR: Failed innerproduct operator*  test." << endl;
		failures++;
	}
	cout << endl;
	
	// testing isColinear
    ntuple<double,3> a3(1.0, 2.0, 3.0);
    ntuple<double,3> b3(2.0, 4.0, 6.0);
    ntuple<double,3> c3(3.0, 6.0, 9.0);
    ntuple<double,3> d3(3.0, 6.0, 9.0 + 1.0E-06);
    ntuple<double,3> f3(2.0, 1.0, 3.0);

    float va6[] = {1.0, 2.0, 3.0, 1.0, 2.0, 0.0};
    float vb6[] = {2.0, 4.0, 6.0, 2.0, 4.0, 0.0};
    float vc6[] = {3.0, 6.0, 9.0, 3.0, 6.0, 0.0};
    float ve6[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0E-02};

    ntuple<float,6> a6(va6);
    ntuple<float,6> b6(vb6);
    ntuple<float,6> c6(vc6);
    ntuple<float,6> e6(ve6);
    ntuple<float,6> d6 = c6 + e6;

    ntuple<int,2> a2(1000, 2000);
    ntuple<int,2> b2(2000, 4000);
    ntuple<int,2> c2(3000, 6000);

    ntuple<int,3> i3(1, 2, 3);
    ntuple<int,3> j3(0, 0, 0);
    ntuple<int,3> k3(-1, -2, -3);

    if (  isColinear(a3, b3, c3, 0.0)    
          && isColinear(a3, a3, a3, 0.0) 
          && ! isColinear(a3, b3, d3, 0.0) 
          && isColinear(a3, b3, d3, 1.0E-03)
          && isColinear(a6, b6, c6, 0.0)    
          && ! isColinear(a6, b6, d6, 0.0) 
          && isColinear(a6, b6, d6, 1.0E-03)
          &&  isColinear(a2, b2, c2, 0.0)
          && isColinear(i3, j3, k3)
       ) {
      cout << "Passed isColinear test." << endl;
      successes++;
    } else {
      failures++;
      cout << "ERROR: Failed isColinear test." << endl;
    }
    cout << endl;

    if (  isColinear3d(a3, b3, c3)    
          && isColinear3d(a3, a3, a3, 0.0) 
          && ! isColinear3d(a3, b3, d3, 0.0) 
          && isColinear3d(a3, b3, d3, 1.0E-03)
          && isColinear3d(i3, j3, k3)
       ) {
      cout << "Passed isColinear3d test." << endl;
      successes++;
    } else {
      failures++;
      cout << "ERROR: Failed isColinear3d test." << endl;
    }
	cout << endl;

    if (  isCoplanar3d(a3, b3, c3, d3, 0.0)
          && ! isCoplanar3d(ihat, jhat, khat, zero3)
          && isCoplanar3d(xhat, yhat, zhat, (xhat + yhat + zhat)/3.0f, 1.0E-07)
          && ! isCoplanar3d(a3, c3, f3, a3 ^ f3)
       ) {
      cout << "Passed isCoplanar3d test." << endl;
      successes++;
    } else {
      failures++;
      cout << "ERROR: Failed isCoplanar3d test." << endl;
    }
	cout << endl;


	// Test modulus.
	ntuple<int,3> a(6,12,18);
	ntuple<int,3> b(1,2,3);
	ntuple<int,3> c(4, 8, 12);

	testStatus = (a.modulus() == b);
	if (! testStatus) {
	  ERROR_LOG("a.modulus() == b");
	  failures++;
	} else {
	  successes++;
	}

	// Validate that modulus is nondestructive
	testStatus = (a == ntuple<int, 3>(6, 12, 18));
	if (! testStatus) {
	  ERROR_LOG("a == ntuple<int, 3>(6, 12, 18)");
	  failures++;
	} else {
	  successes++;
	}

	// Compound test of modulus
	testStatus = (a.modulus() == c.modulus());
	if (! testStatus) {
	  ERROR_LOG("a.modulus() == c.modulus()");
	  failures++;
	} else {
	  successes++;
	}

	ntuple<double, 3> x(1.0, 2.0, 3.0);
	ntuple<double, 3> y(2.0, 4.0, 6.0);
	ntuple<double, 3> xr = x.modulus();
	ntuple<double, 3> yr = y.modulus();
	testStatus = withinEpsilon(xr[0], yr[0]) && withinEpsilon(xr[1], yr[1]) && withinEpsilon(xr[2], yr[2]);
	if (! testStatus) {
	  ERROR_LOG("xr is within epslion of yr");
	  failures++;
	} else {
	  successes++;
	}

	// Test reduce
	a.reduce();
	testStatus = (a == b);
	if (! testStatus) {
	  ERROR_LOG("reduce(a==b)");
	  failures++;
	} else {
	  successes++;
	}

	
	// CODA: Print out the test result summar:
	cout << "Test summary: " << successes << " successes, and " << failures << " failures." << endl;
	return 0;	
}
