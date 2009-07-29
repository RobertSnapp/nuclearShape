/*****************************************************************************
 * voxelImage.h 
 * 
 * Robert R. Snapp  Copyright (C) 2006.
 *****************************************************************************/

#ifndef __VOXEL_IMAGE_H__
#define __VOXEL_IMAGE_H__

#include <algorithm> // min function
#include <cmath> // log function
#include <valarray>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cassert>
#include <typeinfo>
#include "voxelArray.h"
#include "env.h"  // OpenGL declarations
#include "lsm.h"  // for types uint8, etc.

typedef unsigned char uchar;
typedef unsigned char uint8;
typedef unsigned long ulong;
typedef unsigned short ushort;

template<typename T> 
class voxelImageBase : public voxelArray {
protected:
  std::valarray<T> d_voxel;
  uint32 d_bytesPerBand;
  uint32 d_bytesPerLayer;
  uint32 d_bytesPerRow;
	
public:
 voxelImageBase(ulong r, ulong c, ulong l, ulong b) : voxelArray(r, c, l, b) {
	d_bytesPerRow   = d_cols   * sizeof(T);
	d_bytesPerLayer = d_rows   * d_bytesPerRow;
	d_bytesPerBand  = d_layers * d_bytesPerLayer;
	d_voxel.resize(d_rows*d_cols*d_layers*d_bands);		
  }

 voxelImageBase() :  voxelArray(), d_voxel(0), d_bytesPerBand(0), d_bytesPerLayer(0), d_bytesPerRow(0) {}
	
 virtual ~voxelImageBase() {}
		
 voxelImageBase(const voxelImageBase & p) : 
    d_bytesPerBand(p.d_bytesPerBand), d_bytesPerLayer(p.d_bytesPerLayer), d_bytesPerRow(p.d_bytesPerRow), voxelArray(p) {		
	  d_voxel.resize(p.d_rows*p.d_cols*p.d_layers*p.d_bands);
	  d_voxel = p.d_voxel;  // Deep copy.
 }
		
 void clear() {
   voxelArray::clear();
   d_voxel.resize(0);
   d_bytesPerRow = 0;
   d_bytesPerLayer = 0;
   d_bytesPerBand = 0;
 }
	
 voxelImageBase& operator=(const voxelImageBase &p) {
   // guard against self-assignment
   if (this != &p) {
	 static_cast<voxelArray&>(*this) = static_cast<const voxelArray&>(p);
	 d_bytesPerRow   = p.d_bytesPerRow;
	 d_bytesPerLayer = p.d_bytesPerLayer;
	 d_bytesPerBand  = p.d_bytesPerBand;
	 d_voxel.resize(p.d_rows*p.d_cols*p.d_layers*p.d_bands);
	 d_voxel = p.d_voxel;
   }
   return *this;
 }
	
 uint32 bytesPerLayer() const {return d_bytesPerLayer;}
 uint32 bytesPerRow() const {return d_bytesPerRow;}
 uint32 bytesPerBand() const {return d_bytesPerBand;}
	
 T& getVoxel(ulong i, ulong j, ulong l, ulong b) {
   assert(0 <= i && i < d_rows);
   assert(0 <= j && j < d_cols);
   assert(0 <= l && l < d_layers);
   assert(0 <= b && b < d_bands);
   return d_voxel[((b*d_layers + l)*d_rows + i)*d_cols + j];
 }
			
 void setVoxel(ulong i, ulong j, ulong l, ulong b, T v) {
   assert(0 <= i && i < d_rows);
   assert(0 <= j && j < d_cols);
   assert(0 <= l && l < d_layers);
   assert(0 <= b && b < d_bands);
   d_voxel[((b*d_layers + l)*d_rows + i)*d_cols + j] = v;
 }
    
 void getVoxel(ulong i, ulong j, ulong l, ulong b, T &v) const {
   v = d_voxel[((b*d_layers + l)*d_rows + i)*d_cols + j];
 }
    
 T getVoxelClip(ulong i, ulong j, ulong l, ulong b) const {
   if (0 <= i && i < d_rows && 
	   0 <= j && j < d_cols && 
	   0 <= l && l < d_layers &&
	   0 <= b && b < d_bands    ) {
	 return d_voxel[((b*d_layers + l)*d_rows + i)*d_cols + j];
   } else {
	 return 0;
   }
 }
    
 T* getVoxelOffset(ulong i, ulong j, ulong l, ulong b) {
   return &d_voxel[((b*d_layers + l)*d_rows + i)*d_cols + j];
 }

 void resize() {
   d_bytesPerRow   = d_cols   * sizeof(T);
   d_bytesPerLayer = d_rows   * d_bytesPerRow;
   d_bytesPerBand  = d_layers * d_bytesPerLayer;
   d_voxel.resize(d_rows*d_cols*d_layers*d_bands);	
 }
	
 void setColorMap(bool inverseVideo) {
   	if (inverseVideo) {
	  glPixelTransferf(GL_RED_SCALE,   -1.0);
	  glPixelTransferf(GL_GREEN_SCALE, -1.0);
	  glPixelTransferf(GL_BLUE_SCALE,  -1.0);
	  glPixelTransferf(GL_RED_BIAS,     1.0);
	  glPixelTransferf(GL_GREEN_BIAS,   1.0);
	  glPixelTransferf(GL_BLUE_BIAS,    1.0);
	} else {
	  glPixelTransferf(GL_RED_SCALE,    1.0);
	  glPixelTransferf(GL_GREEN_SCALE,  1.0);
	  glPixelTransferf(GL_BLUE_SCALE,   1.0);
	  glPixelTransferf(GL_RED_BIAS,     0.0);
	  glPixelTransferf(GL_GREEN_BIAS,   0.0);
	  glPixelTransferf(GL_BLUE_BIAS,    0.0);
	}
	return;
 }

  // function drawPixels renders a rectangular subregion of the indicated layer and band. The arguments
  // width and height specify the dimensions of the subregion, and col and row specify the location
  // of the lower left-hand corner of the subregion within the full layer. Note that the defaults are
  // defined such that the call drawPixels(layer, band) should display the entire layer.
  virtual void drawPixels(uint32 layer, uint32 band, 
                          uint32 width = 0, uint32 height = 0, 
                          uint32 col = 0, uint32 row = 0, bool inverseVideo = false) = 0; 
  virtual void drawPixels(int layer, int band, 
                          int width = 0, int height = 0, 
                          int col = 0, int row = 0, bool inverseVideo = false) = 0; 
};


template<typename T>
class voxelImage : public voxelImageBase<T> {
public:
	voxelImage<T>() : voxelImageBase<T>() {}
	voxelImage<T>(ulong r, ulong c, ulong l, ulong b) : voxelImageBase<T>(r, c, l, b) {}
};

template<>
class voxelImage<uint8> : public voxelImageBase<uint8> {
public:
  voxelImage<uint8>() : voxelImageBase<uint8>() {}
  voxelImage<uint8>(ulong r, ulong c, ulong l, ulong b) : voxelImageBase<uint8>(r, c, l, b) {}
  void drawPixels(uint32 layer, uint32 band, 
				  uint32 width = 0, uint32 height = 0, 
				  uint32 col0 = 0, uint32 row0 = 0,
				  bool inverseVideo = false) {
    assert(0 <= row0 && row0 < d_rows);
    assert(0 <= col0 && col0 < d_cols);
    assert(0 <= width && width <= d_cols - col0);
    assert(0 <= height && height <= d_rows - row0);

    if (width == 0) {
      width = d_cols - col0;
    }
    if (height == 0) {
      height = d_rows - row0;
    }

	setColorMap(inverseVideo);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, d_cols);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, row0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, col0);
    glDrawPixels(static_cast<GLsizei>(width), 
                 static_cast<GLsizei>(height), 
                 GL_LUMINANCE, GL_UNSIGNED_BYTE, 
                 &d_voxel[(band*d_layers + layer)*d_bytesPerLayer]);

    // Reset to default values
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  }

  void drawPixels(int layer, int band, int width = 0, int height = 0, int col0 = 0, int row0 = 0, bool inverseVideo = false) {
     drawPixels(static_cast<uint32>(layer), 
                static_cast<uint32>(band),
                static_cast<uint32>(width),
                static_cast<uint32>(height),
                static_cast<uint32>(col0),
                static_cast<uint32>(row0),
				inverseVideo);
   }
};

#ifdef COMMENT
template<>
void voxelImage<GLshort>::drawPixels<GLshort>(uint32 layer, uint32 band) {
	glDrawPixels(static_cast<GLsizei>(d_cols), 
				 static_cast<GLsizei>(d_rows), 
				 GL_LUMINANCE, GL_SHORT, 
				 &d_voxel[(band*d_layers + layer)*d_bytesPerLayer]);
}

template<>
void voxelImage<GLint>::drawPixels<GLint>(uint32 layer, uint32 band) {
	glDrawPixels(static_cast<GLsizei>(d_cols), 
				 static_cast<GLsizei>(d_rows), 
				 GL_LUMINANCE, GL_INT, 
				 &d_voxel[(band*d_layers + layer)*d_bytesPerLayer]);
}

template<>
void voxelImage<GLbyte>::drawPixels<GLbyte>(uint32 layer, uint32 band) {
	glDrawPixels(static_cast<GLsizei>(d_cols), 
				 static_cast<GLsizei>(d_rows), 
				 GL_LUMINANCE, GL_BYTE, 
				 &d_voxel[(band*d_layers + layer)*d_bytesPerLayer]);
}
#endif

#ifdef COMMENT
template<>
class grayImage<GLshort> : public grayBase<GLshort> {
public:
	grayImage<GLshort>() : grayBase<GLshort>() {}
	grayImage<GLshort>(ulong r, ulong c) : grayBase<GLshort>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_SHORT, &pixel_[0]);
	}
};

template<>
class grayImage<GLint> : public grayBase<GLint> {
public:
	grayImage<GLint>() : grayBase<GLint>() {}
	grayImage<GLint>(ulong r, ulong c) : grayBase<GLint>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_INT, &pixel_[0]);
	}
};

template<>
class grayImage<GLfloat> : public grayBase<GLfloat> {
public:
	grayImage<GLfloat>() : grayBase<GLfloat>() {}
	grayImage<GLfloat>(ulong r, ulong c) : grayBase<GLfloat>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_FLOAT, &pixel_[0]);
	}
};

template<>
class grayImage<GLdouble> : public grayBase<GLdouble> {
public:
	grayImage<GLdouble>() : grayBase<GLdouble>() {}
	grayImage<GLdouble>(ulong r, ulong c) : grayBase<GLdouble>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_DOUBLE, &pixel_[0]);
	}
};

template<>
class grayImage<GLbyte> : public grayBase<GLbyte> {
public:
	//grayImage<GLbyte>() : grayBase<GLbyte>() {}
	//grayImage<GLbyte>(ulong r, ulong c) : grayBase<GLbyte>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_BYTE, &pixel_[0]);
	}
};

template<>
class grayImage<GLuint> : public grayBase<GLuint> {
public:
	//grayImage<GLuint>() : grayBase<GLuint>() {}
	//grayImage<GLuint>(ulong r, ulong c) : grayBase<GLuint>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_INT, &pixel_[0]);
	}
};

template<>
class grayImage<GLushort> : public grayBase<GLushort> {
public:
	//grayImage<GLushort>() : grayBase<GLushort>() {}
	//grayImage<GLushort>(ulong r, ulong c) : grayBase<GLushort>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_SHORT, &pixel_[0]);
	}
};

#endif

#endif /* __VOXEL_IMAGE__ */
