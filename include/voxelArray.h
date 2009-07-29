/*****************************************************************************
 * voxelArray.h
 * 
 * voxelArray represents a base class for supporting the storage, management, 
 * and manipulation 3d raster images.  This class implements only the most common
 * functionality: management of the image geometry, the title, and possibly
 * other ancillary information (author, creator, filename, etc.)
 *
 * Robert R. Snapp Copyright (C) 2006.
 *****************************************************************************/

#ifndef __VOXELARRAY_H__
#define __VOXELARRAY_H__

#include <string>
#include <cassert>
#include <tiff.h>

// typedef unsigned int uint32;

class voxelArray {
protected:
	uint32 d_rows;
	uint32 d_cols;
	uint32 d_layers;
	uint32 d_bands;

	std::string d_title;
	
	double d_dx; // physical displacement between adjacent columns.
	double d_dy; // physical displacement between adjacent rows.
	double d_dz; // physical displacement between adjacent layers.

public:
 voxelArray() : 
  d_rows(0), d_cols(0), d_layers(0), d_bands(0), d_dx(0), d_dy(0), d_dz(0) {}

 voxelArray(uint32 r, uint32 c, uint32 l, uint32 d) : 
  d_rows(r), d_cols(c), d_layers(l), d_bands(d), d_dx(0), d_dy(0), d_dz(0) {
	assert(r > 0);
	assert(c > 0);
	assert(l > 0);
	assert(d > 0);
  };
 virtual ~voxelArray() {};


 voxelArray(const voxelArray &v) : 
  d_rows(v.d_rows), d_cols(v.d_cols), d_layers(v.d_layers), d_bands(v.d_bands), d_title(v.d_title),
	d_dx(v.d_dx), d_dy(v.d_dy), d_dz(v.d_dz) { };


  voxelArray& operator=(const voxelArray &v) {
	// guard against self-assignment
	if (this != &v) {
	  d_rows   = v.d_rows;
	  d_cols   = v.d_cols;
	  d_layers = v.d_layers;
	  d_bands  = v.d_bands;
	  d_title  = v.d_title;
	  d_dx = v.d_dx;
	  d_dy = v.d_dy;
	  d_dz = v.d_dz;

	}	
	return *this;
  }
	
   	uint32 rows() const {return d_rows;}
   	uint32 cols() const {return d_cols;}
   	uint32 layers() const {return d_layers;}
   	uint32 bands() const {return d_bands;}

    std::string title() const {
    	return d_title;
    }
    
   void setRows(int r) {
    	assert(r > 0); 
    	d_rows = r;
    }
    void setCols(int c) {
    	assert(c > 0); 
    	d_cols = c;
    }

	void setLayers(int l) {
    	assert(l > 0); 
    	d_layers = l;
    }

	void setBands(int b) {
    	assert(b > 0); 
    	d_bands = b;
    }

    void setTitle(std::string t) {
    	d_title = t;
    }
    
    void clear() {
    	d_rows = 0;
    	d_cols = 0;
    	d_bands = 0;
		d_layers = 0;
    	d_title.clear();
		d_dx     = 0;
		d_dy     = 0;
		d_dz     = 0;
    }
    
	double getDX() const {return d_dx;}
	double getDY() const {return d_dy;}
	double getDZ() const {return d_dz;}
	
	void setDX(double dx) {d_dx = dx;}
	void setDY(double dy) {d_dy = dy;}
	void setDZ(double dz) {d_dz = dz;}
	
};


#endif /* __VOXELARRAY_H__ */
