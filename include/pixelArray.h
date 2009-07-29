/*****************************************************************************
 * pixelArray.h
 * 
 * pixelArray represents a base class for supporting the storage, management, 
 * and manipulation raster images.  This class implements only the most common
 * functionality: management of the image geometry, the title, and possibly
 * other ancillary information (author, creator, filename, etc.)
 *
 * Robert R. Snapp Copyright (C) 2006.
 *****************************************************************************/

#ifndef __PIXELARRAY_H__
#define __PIXELARRAY_H__

#include <iostream>
#include <ostream>
#include <string>
#include <cassert>
#include <tiff.h>  // defines the type uint32

// typedef unsigned long int uint32;

class pixelArray {
 protected:
  uint32 d_rows;
  uint32 d_cols;
  uint32 d_bands;
  uint32 d_bytesPerRow;
  std::string d_title;

public:
 pixelArray() : d_rows(0), d_cols(0), d_bands(0), d_bytesPerRow(0) {}
 pixelArray(uint32 r, uint32 c, uint32 d) : d_rows(r), d_cols(c), d_bands(d), d_bytesPerRow(0) {};

 virtual ~pixelArray() {};

 // Assume default copy constructor

 // Assignment operator
  pixelArray& operator=(const pixelArray &p) {
    if (this != &p) {
      d_rows  = p.d_rows;
      d_cols  = p.d_cols;
      d_bands = p.d_bands;
	  d_bytesPerRow = p.d_bytesPerRow;
      d_title = p.d_title;
    }
    return *this;
  }
	
  uint32 rows() {return d_rows;}
  uint32 cols() {return d_cols;}
  uint32 bands() {return d_bands;}
  uint32 bytesPerRow() {return d_bytesPerRow;}
    
  std::string title() const {
    return d_title;
  }
    
  void setTitle(std::string s) {
    d_title = s;
  }
    
  void clear() {
    d_rows = 0;
    d_cols = 0;
    d_bands = 0;
    d_bytesPerRow = 0;
    d_title.clear();
  }
  
  void describe(std::ostream &os) const {
    os << "d_rows = " << d_rows << ", "
       << "d_cols = " << d_cols << ", "
       << "d_bands = " << d_bands << ", "
       << "d_bytesPerRow = " << d_bytesPerRow << ", "
       << "d_title = " << d_title
       << std::endl;
  }
};


#endif /* __PIXELARRAY_H__ */
