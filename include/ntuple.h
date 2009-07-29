/*** ntuple.h ***/

#ifndef __NTUPLE__H__
#define __NTUPLE__H__

/* Implements an n-dimensional vector class for computer graphics applications. It is called
 * ntuple to distinguish it clearly from the standard template vector, and valarray.
 *
 * Robert R. Snapp
 * Copyright 2006
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>


// forward declarations:
template<typename T, size_t N> class ntuple;

template<typename T, size_t N>
ntuple<T,N> operator*(const T scalar, const ntuple<T,N>& v);

template<typename T, size_t N> 
ntuple<T,N> operator*(const ntuple<T,N>& v, const T scalar);

template<typename T, size_t N> 
ntuple<T,N> operator/(const ntuple<T,N>& v, const T scalar);
	
template<typename T, size_t N> 
class ntuple {
 private:
  T d_s[N];  // store

 public:
  typedef T value_type;
  static bool const debug_ntuple_h = false; // This flag can be set to true for debugging.
	
	// An explicit constructor that allows specifying up to at most the first four components.
	explicit ntuple<T,N>(const T& xp = T(0), 
				const T& yp = T(0), 
				const T& zp = T(0),
				const T& wp = T(0)) {
		if (N > 0) {
			d_s[0] = xp;
			if (N > 1) {
			  d_s[1] = yp;
			  if (N > 2) {
				d_s[2] = zp;
				if (N > 3) {
				  d_s[3] = wp;
				  for (size_t i = 4; i < N; i++)
					d_s[i] = T(0);	// Any additional components are initialized to zero.
				}}}}
	}
	
	// A constructor that initializes values from a standard array.
	explicit ntuple<T,N>(const T* v) {
		for (size_t i = 0; i < N; i++) {
			d_s[i] = v[i];
		}
	}
	
	// Copy constructor.
	ntuple<T,N>(const ntuple<T,N>& v) {
	  
		for (size_t i = 0; i < N; i++)
			d_s[i] = v.d_s[i];
	}
	
	// Assign the value of an ntuple from another ntuple of the same type.
	ntuple<T,N>& operator=(const ntuple<T,N>& v) {
	  // guard against self-assignment
      if (this != &v) {
		for (size_t i = 0; i < N; i++) {
		  d_s[i] = v.d_s[i];
		}

      }
      return *this;
	}
	
	~ntuple<T,N>() {}
			
	ntuple<T,N>& clear() {
		for (register size_t i = 0; i < N; i++) {
			d_s[i] = 0;
		}
		return *this;
	}
	
	bool operator==(const ntuple<T,N>& v1) const {
		for(size_t i = 0; i < N; i++) {
			if (d_s[i] != v1[i]) return false;
		}
		return true;
	}
	
	bool operator!=(const ntuple<T,N>& v1) const {
		for(size_t i = 0; i < N; i++) {
			if (d_s[i] != v1[i]) return true;
		}
		return false;
	}
		
	bool operator<(const ntuple<T,N>& v1) const {
		for (register size_t i = 0; i < N; i++) {
			if(! (d_s[i] < v1[i])) return false;
		}
		return true;
	}	

	ntuple<T,N> operator+(const ntuple<T,N>& v1) const {
		ntuple<T,N> sum;
		for (size_t i = 0; i < N; i++)
			sum[i] = d_s[i] + v1[i];
		return sum;
	}
	
	// friend 
	// ntuple<T,N> operator-(const ntuple<T,N>& v1, const ntuple<T,N>& v2) {
	// 	ntuple<T,N> diff;
	// 	for (size_t i = 0; i < N; i++)
	// 		diff.d_s[i] = v1.d_s[i] - v2.d_s[i];
	// 	return diff;
	// }
	
	ntuple<T,N> operator-(const ntuple<T,N>& v1) const {
		ntuple<T,N> diff;
		for (size_t i = 0; i < N; i++)
			diff[i] = d_s[i] - v1[i];
		return diff;
	}
	
	friend ntuple<T,N> operator*<>(const T scalar, const ntuple<T,N>& v);	
	friend ntuple<T,N> operator*<>(const ntuple<T,N>& v, const T scalar);
	friend ntuple<T,N> operator/<>(const ntuple<T,N>& v, const T scalar);

	
	ntuple<T,N>& operator+=(const ntuple<T,N>& v1) {
		for (size_t i = 0; i < N; i++) {
			d_s[i] += v1.d_s[i];
		}
		return *this;
	}
	
	ntuple<T,N>& operator-=(const ntuple<T,N>& v1) {
		for (size_t i = 0; i < N; i++) {
			d_s[i] -= v1.d_s[i];
		}
		return *this;
	}
	
	ntuple<T,N>& operator*=(const T scalar) {
		for (size_t i = 0; i < N; i++) {
			d_s[i] *= scalar;
		}
		return *this;
	}
	
	ntuple<T,N>& operator/=(const T scalar) {
		for(size_t i = 0; i < N; i++) {
			d_s[i] /= scalar;
		}
		return *this;
	}
	
	template<typename T1, size_t N1>
	double innerProduct(const ntuple<T1,N1> & v) const {
		double value = 0;
		for(size_t i = 0; i < std::min(N,N1); i++) {
			value += d_s[i]*v[i];
		}
		return value;
	}
	
	ntuple<T,N> interpolate(const ntuple<T,N>&v, double f = 0.5) const {
		ntuple<T,N> val;
		val = (1-f)*(*this) + f*v;
		return val;
	}

#ifdef COMMENT
	template<typename T1, typename T2>
	friend double innerProduct (const ntuple<T1,N>& v, const ntuple<T2,N>& wp) {
		double value = double(0);
		for (size_t i = 0; i < N; i++) {
			value += v.d_s[i] * wp.d_s[i];
		}
		return value;
	}
	
	template<typename T1, typename T2>
	friend ntuple<T,N> interpolate(const ntuple<T1,N>&v, const ntuple<T2,N>&w, double f) {
		ntuple<T,N> val;
		val = (1-f)*v + f*w;
		return val;
	}
#endif
	
	void swap(ntuple<T,N>& v) {
		for(size_t i = 0; i < N; i++) {
			T temp = d_s[i];
			d_s[i] = v.d_s[i];
			v.d_s[i] = temp;
		}
	}
	
	double squareNorm(void) const {
		return  (*this) * (*this);  // inner product
	}
	
	double norm(void) const {
		return sqrt(this->squareNorm());
	}
	
	double infinityNorm() {
		double maxComponent = 0;
		for(size_t i = 0; i < N; i++) {
			maxComponent = std::max(maxComponent, std::abs(d_s[i]));
		}
		return maxComponent;
	}
	
	double pNorm(double p = 2.0) {
		double sum = 0.0;
		for(size_t i = 0; i < N; i++) {
			sum += std::pow(static_cast<double>(std::abs(d_s[i])), p);
		}
		return std::pow(sum, 1.0/p);
	}
	
	ntuple<T,N> normalize(double p = 2.0) {
		return *this /= pNorm(p);
	}
	
	// Reduce is a generalization of normalize, in that it rescales the length
	// of each ntuple, by dividing each component by the generalized greatest
	// common divisor of 
	ntuple<T,N> reduce() {
	  return *this /= genGCD(*this);
	}

	bool isWithinEpsilon(const ntuple<T,N>& center, double epsilon, double p = 2.0) const {
		ntuple<T,N> diff = *this - center;
		if (p == 2.0) {
			return diff.squareNorm() <= epsilon*epsilon;
		} else {
			return diff.pNorm(p) <= epsilon;
		}
	}
		
	T operator[](size_t i) const {
		assert(0 <= i && i < N);
		return d_s[i];
	}
	
	T& operator[](size_t i) {
		assert(0 <= i && i < N);
		return d_s[i];
	}

	T* array() {
		return &d_s[0];
	}
	
	// For double and float types, modulus returns a normalized n-tuple. For other types, each component
	// is reduced by the greatest common divisor of the n-tuple coordinates.
	ntuple<T,N> modulus() const {
	  T d = genGCD(*this);
	  return *this/d;
	}
	  
	T x() const {return d_s[0];}
	T y() const {return (N > 1 ? d_s[1] : T());}
	T z() const {return (N > 2 ? d_s[2] : T());}
	T w() const {return (N > 3 ? d_s[3] : T());}
	
	T& x()  {return d_s[0];}
	T& y()  {return (N > 1 ? d_s[1] : T());}
	T& z()  {return (N > 2 ? d_s[2] : T());}
	T& w()  {return (N > 3 ? d_s[3] : T());}
	
	
	void x(const T& xp) {d_s[0] = xp;}
	void y(const T& yp) {if (N > 1) d_s[1] = yp;}
	void z(const T& zp) {if (N > 2) d_s[2] = zp;}
	void w(const T& wp) {if (N > 3) d_s[3] = wp;}
		

};

template<typename T>
T gcd(T u, T v) {
  T r;
  u = std::abs(u);
  v = std::abs(v);
  while (v > 0) {
	r = u % v;
	u = v;
	v = r;
  }
  return u;
}

template<size_t N>
double genGCD(ntuple<double, N> v) {
  return v.norm();
}

template<size_t N>
float genGCD(ntuple<float, N> v) {
  return v.norm();
}


template<typename T, size_t N>
T genGCD(ntuple<T,N> v) {
  T      d = v[0];
  size_t k = 1;

  while (d != 1 && k < N) {
	d = gcd<T>(d, v[k]);
	k++;
  }

#ifdef COMMENT  
  T 	 d = v[N-1];
  size_t k = N - 1;

  while (d != 1 && k >= 0) {
	d = gcd<T>(v[k], d);
	k--;
  }
#endif

  return d;
}
  
// scalar-vector product
template<typename T, size_t N>
ntuple<T,N> operator*(const T scalar, const ntuple<T,N>& v) {
	ntuple<T,N> prod;
	for (size_t i = 0; i < N; i++)
		prod.d_s[i] = scalar*v.d_s[i];
	return prod;
}

// vector-scalar product
template<typename T, size_t N>
ntuple<T,N> operator*(const ntuple<T,N>& v, const T scalar) {
	ntuple<T,N> prod;
	for (size_t i = 0; i < N; i++)
		prod.d_s[i] = scalar*v.d_s[i];
	return prod;
}

// inner-product
template<typename T, size_t N1, size_t N2>
T operator*(const ntuple<T,N1>& xp, const ntuple<T,N2>& yp) {
	T value = 0;
	for(size_t i = 0; i < std::min(N1,N2); i++) {
		value += xp[i]*yp[i];
	}
	return value;
}

#ifdef COMMENT
template<typename T1, typename T2, size_t N>
  bool operator==(const ntuple<T1,N> &x, const ntuple<T2,N> &y) {
  for (size_t i = 0; i < N; i++) {
	if (x[i] != y[i]) return false;
  }
  return true;
}
#endif

template<typename T1, typename T2, size_t N1, size_t N2>
double operator*(const ntuple<T1,N1>& xp, const ntuple<T2,N2>& yp) {
	double value = 0;
	for(size_t i = 0; i < std::min(N1,N2); i++) {
		value += xp[i]*yp[i];
	}
	return value;
}


template<typename T, size_t N>
ntuple<T,N> operator/(const ntuple<T,N>& v, const T scalar) {
	ntuple<T,N> quotient;
	for (size_t i = 0; i < N; i++)
		quotient.d_s[i] = v.d_s[i]/scalar;
	return quotient;
}

template<typename T>
inline ntuple<T,2> perp(ntuple<T,2> const &xp) {
  ntuple<T,2> result;
  result[0] =   xp[1];
  result[1] = - xp[0];
  return result;
}

template<typename T>
ntuple<T,3> crossProduct(const ntuple<T,3> &xp, const ntuple<T,3> &yp) {
	ntuple<T,3> result;
	result[0] = xp[1]*yp[2] - xp[2]*yp[1];
	result[1] = xp[2]*yp[0] - xp[0]*yp[2];
	result[2] = xp[0]*yp[1] - xp[1]*yp[0];
	return result;
}

// Computes the three-dimensional cross product, for operands with common type.
template<typename T>
ntuple<T, 3> operator^(const ntuple<T,3> &xp, const ntuple<T,3> &yp) {
	ntuple<T,3> result;
	result[0] = xp[1]*yp[2] - xp[2]*yp[1];
	result[1] = xp[2]*yp[0] - xp[0]*yp[2];
	result[2] = xp[0]*yp[1] - xp[1]*yp[0];
	return result;
}

template<typename T1, typename T2>
ntuple<double, 3> operator^(const ntuple<T1,3> &xp, const ntuple<T2,3> &yp) {
	ntuple<double,3> result;
	result[0] = xp[1]*yp[2] - xp[2]*yp[1];
	result[1] = xp[2]*yp[0] - xp[0]*yp[2];
	result[2] = xp[0]*yp[1] - xp[1]*yp[0];
	return result;
}

template<typename T>
T tripleProduct(const ntuple<T,3> &a, const ntuple<T,3> &b, const ntuple<T,3> &c) {
	return  a[0]*(b[1]*c[2] - b[2]*c[1]) + 
			a[1]*(b[2]*c[0] - b[0]*c[2]) +
			a[2]*(b[0]*c[1] - b[1]*c[0]);
}

// Function isColinear3d tests if three, 3-dimensional ntuples are colinear by computing the two-dimensional area
// of the triangular region they delimit. If 4*area^2 does not exceed the product of eps multiplied by a geometric
// scale factor, then the predicate returns true; otherwise false.
template<typename T>
bool isColinear3d(const ntuple<T,3> &a, const ntuple<T,3> &b, const ntuple<T,3> &c, double eps = 0) {
  ntuple<T,3> d1 = b - a;                                   // side 1
  ntuple<T,3> d2 = c - a;                                   // side 2
  ntuple<T,3> cp = d1 ^ d2;                                 // a vector whoops, only defined for three 3.
  T ip = cp * cp;                                           // 4 * area^2
  double sf = static_cast<double>(max(d1*d1, d2*d2));       // scale factor

  if (ntuple<T,3>::debug_ntuple_h) {
    std::cout << "a = " << std::setprecision(15) << a << ", b = " << b << ", c = " << c 
              << ", lhs = " << abs(ip) 
              << ", rhs = " <<  eps*sf << std::endl;
  }

  return static_cast<double>(abs(ip)) <= eps*sf;
}

// Function isCollinear tests if three, n-dimensional ntuples are colinear by computing the angle between
// two of the rays defined by these three points. Consequently, only inner products are required.
template<typename T, size_t N>
inline bool isColinear(const ntuple<T,N> &a, const ntuple<T,N> &b, const ntuple<T,N> &c, double eps = 0) {
  ntuple<T,N> d1 = b - a;                                   // side 1
  ntuple<T,N> d2 = c - a;                                   // side 2
  T ip  = d1 * d2;
  T ip2 = ip * ip;
  T sf = (d1*d1) * (d2*d2);

  if (ntuple<T,N>::debug_ntuple_h) {
    std::cout << "a = " << std::setprecision(15) << a << ", b = " << b << ", c = " << c 
              << ", lhs = " << (1.0 - eps)*static_cast<double>(sf) 
              << ", rhs = " <<  static_cast<double>(ip2) << std::endl;
  }

  return sf == 0 || (1.0 - eps)*static_cast<double>(sf) <= static_cast<double>(ip2);
}

template<typename T>
inline bool isCoplanar3d(const ntuple<T,3> &a, const ntuple<T,3> &b, 
						 const ntuple<T,3> &c, const ntuple<T,3> &d, double eps = 0) {
  ntuple<T,3> d1 = b - a;
  ntuple<T,3> d2 = c - a;
  ntuple<T,3> d3 = d - a;
  T tp = tripleProduct(d1, d2, d3);

  if (ntuple<T,3>::debug_ntuple_h) {
    std::cout << "a = " << std::setprecision(15) << a 
			  << ", b = " << b 
			  << ", c = " << c 
			  << ", d = " << d
              << ", lhs = " << abs(tp) 
              << ", rhs = " <<  eps << std::endl;
  }


  return static_cast<double>(abs(tp)) <= eps;
}
  
// Read in a ntuple object from a stream that is delineated by "(", ",", and ")", such
// as (1.25, 3.24), or (32, -45).
template<typename T, size_t N, typename charT, typename traits>
std::basic_istream<charT, traits>& 
operator>>(std::basic_istream<charT, traits>& is, ntuple<T,N>& v) {
	T buffer;
	char c;
	
	if (! (is >> c)) return is;
	if (c != '(') {
		is.setstate(std::ios_base::failbit);
		return is;
	}
	
	for (size_t i = 0; i < N - 1; i++) {
		if (! (is >> buffer)) return is;
		v[i] = buffer;
		
		if (! (is >> c)) return is;
		if (c != ',') {
			is.setstate(std::ios_base::failbit);
			return is;
		}
	}
	
	if (! (is >> buffer)) return is;
		v[N-1] = buffer;
	
	
	if (! (is >> c)) return is;
	if (c != ')') {
		is.setstate(std::ios_base::failbit);
		return is;
	}
	
	return is;
}
	
template<typename T, size_t N, typename charT, typename traits>
std::basic_ostream<charT, traits>& 
operator<<(std::basic_ostream<charT, traits>& os, const ntuple<T,N>& v) {
	os << "(";
	
	for (size_t i = 0; i < N - 1; i++) 
		os << v[i] << ", ";
		
	os << v[N-1] << ")";
	return os;
}


#endif
