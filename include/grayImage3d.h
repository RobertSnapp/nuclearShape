/*****************************************************************************/
/* grayImage3d.h 
/* 
/* Robert R. Snapp  Copyright (C) 2006.
/*****************************************************************************/

#ifndef __GRAYIMAGE3D_H__
#define __GRAYIMAGE3D_H__

#include <algorithm> // min function
#include <cmath> // log function
#include <valarray>
#include <GLUT/glut.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cassert>
#include <typeinfo>
#include "voxelArray.h"


typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;

template<typename T> 
class grayBase3d : public voxelArray {
protected:
	std::valarray<T> d_voxel;
	int d_bytesPerRow;
	int d_bytesPerLayer;
public:
	grayBase3d(int rows, int cols, int layers) : voxelArray(rows, cols, layers, 1) {
		d_bytesPerRow = d_cols*sizeof(T);
		d_bytesPerLayer = d_rows*d_cols*sizeof(T);	
		d_voxel.resize(d_rows*d_cols*d_layers);		
	}

 grayBase3d() : voxelArray(), d_bytesPerRow(0), d_bytesPerLayer(0) {
		voxelArray::d_bands = 1;
	}
	
	~grayBase3d() {}
		
	grayBase3d(const grayBase3d & p) : voxelArray(p) {		
		d_voxel.resize(d_rows*d_cols*d_layers);
		d_voxel = p.d_voxel;  // Deep copy.
		d_bytesPerRow   = p.d_bytesPerRow;
		d_bytesPerLayer = p.d_bytesPerLayer;
	}
		
	void clear() {
		voxelArray::clear();
		d_bytesPerRow   = 0;
		d_bytesPerLayer = 0;
		d_voxel.resize(0);
	}
	
	grayBase3d& operator=(const grayBase3d &p) {
		// guard against self-assignment
		if (this != &p) {
			static_cast<voxelArray&>(*this) = static_cast<const voxelArray&>(p);
			d_voxel.resize(d_rows*d_cols*d_layers);
			d_voxel = p.d_voxel;
		}
		return *this;
	}
	
	T& getVoxel(ulong i, ulong j, ulong l) {
		assert(0 <= i); assert(i < d_rows);
		assert(0 <= j); assert(j < d_cols);
		assert(0 <= l); assert(l < d_layers);
		return d_voxel[l*d_rows*d_cols + i*d_cols + j];
	}
			
    void setVoxel(ulong i, ulong j, ulong l, T v) {
    	assert(0 <= i); assert(i < d_rows);
		assert(0 <= j); assert(j < d_cols);
		assert(0 <= l); assert(l < d_layers);
    	d_voxel[l*d_rows*d_cols + i*d_cols + j] = v;
    }
    
    void getVoxel(ulong i, ulong j, ulong l, T &v) const {
    	v = d_voxel[l*d_rows*d_cols + i*d_cols + j];
    }
    
    T getVoxelClip(ulong i, ulong j, ulong l) const {
	  T value = (0 <= i && i < d_rows && 0 <= j && j < d_cols && 0 <= l && l < d_layers ?
				 d_voxel[(l*d_rows + i)*d_cols + j] : 0);
	  return value;
    }
    
 	T* getVoxelOffset(ulong i, ulong j, ulong l) {
		return &d_voxel[l*d_rows*d_cols+ i*d_cols + j];
	}

   
    void resize() {
    	d_bands = 1;
    	d_voxel.resize(d_rows*d_cols*d_layers);
		d_bytesPerRow = d_cols*sizeof(T);
		d_bytesPerLayer = d_rows*d_bytesPerRow;
	}

    void setColorMap(bool inverseVideo) {
	  // glPixelTransferb(GL_MAP_COLOR, true);
	  if (inverseVideo) {
		glPixelTransferf(GL_RED_SCALE,   -1.0);
		glPixelTransferf(GL_GREEN_SCALE, -1.0);
		glPixelTransferf(GL_BLUE_SCALE, -1.0);
		glPixelTransferf(GL_RED_BIAS,   1.0);
		glPixelTransferf(GL_GREEN_BIAS, 1.0);
		glPixelTransferf(GL_BLUE_BIAS,  1.0);
	  } else {
		glPixelTransferf(GL_RED_SCALE,   1.0);
		glPixelTransferf(GL_GREEN_SCALE, 1.0);
		glPixelTransferf(GL_BLUE_SCALE,  1.0);
		glPixelTransferf(GL_RED_BIAS,    0.0);
		glPixelTransferf(GL_GREEN_BIAS,  0.0);
		glPixelTransferf(GL_BLUE_BIAS,   0.0);
	  }
	  return;
	}

};


template<typename T>
class grayImage3d : public grayBase3d<T> {
 public:
  grayImage3d<T>() : grayBase3d<T>() {}
  grayImage3d<T>(ulong r, ulong c, ulong l) : grayBase3d<T>(r, c, l) {}
};

template<>
class grayImage3d<GLubyte> : public grayBase3d<GLubyte> {
public:
  grayImage3d<GLubyte>() : grayBase3d<GLubyte>() {}
  grayImage3d<GLubyte>(const grayImage3d<GLubyte> &i) : grayBase3d<GLubyte>(static_cast<grayBase3d<GLubyte> >(i)) {}
  grayImage3d<GLubyte>(ulong r, ulong c, ulong l) : grayBase3d<GLubyte>(r, c, l) {}

  grayImage3d<GLubyte>& operator=(const grayImage3d<GLubyte> &i) {
    if (this != &i) {
      grayBase3d<GLubyte>::operator=(static_cast<grayBase3d<GLubyte> >(i));
    }
    return *this;
  }

  void drawPixels(ulong l, bool inverseVideo = false) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glDrawPixels(static_cast<GLsizei>(d_cols), 
                 static_cast<GLsizei>(d_rows), GL_LUMINANCE, GL_UNSIGNED_BYTE, &d_voxel[l*d_rows*d_cols]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }

  void rgbPixels(ulong l) {
	assert(l < d_layers);
	setColorMap(false);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glDrawPixels(static_cast<GLsizei>(d_cols), 
                 static_cast<GLsizei>(d_rows), GL_RGB, GL_UNSIGNED_BYTE_3_3_2, &d_voxel[l*d_rows*d_cols]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }
};

template<>
class grayImage3d<GLshort> : public grayBase3d<GLshort> {
public:
  grayImage3d<GLshort>() : grayBase3d<GLshort>() {}
  grayImage3d<GLshort>(const grayImage3d<GLshort> &i) : grayBase3d<GLshort>(static_cast<grayBase3d<GLshort> >(i)) {}
  grayImage3d<GLshort>(ulong r, ulong c, ulong l) : grayBase3d<GLshort>(r, c, l) {}

  grayImage3d<GLshort>& operator=(const grayImage3d<GLshort> &i) {
    if (this != &i) {
      grayBase3d<GLshort>::operator=(static_cast<grayBase3d<GLshort> >(i));
    }
    return *this;
  }
  void drawPixels(ulong l, bool inverseVideo = false) {
	assert(l < d_layers);
	setColorMap(inverseVideo);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_SHORT, &d_voxel[l*d_rows*d_cols]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }
};

template<>
class grayImage3d<GLint> : public grayBase3d<GLint> {
public:
  grayImage3d<GLint>() : grayBase3d<GLint>() {}
  grayImage3d<GLint>(ulong r, ulong c, ulong l) : grayBase3d<GLint>(r, c, l) {}
  void drawPixels(ulong l, bool inverseVideo = false) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_INT, &d_voxel[l*d_rows*d_cols]);
  }

  void mapPixels(ulong l) {
	assert(l < d_layers);
    glDrawPixels(static_cast<GLsizei>(d_cols), 
                 static_cast<GLsizei>(d_rows), GL_COLOR_INDEX, GL_INT, &d_voxel[l*d_rows*d_cols]);
  }
};

template<>
class grayImage3d<GLfloat> : public grayBase3d<GLfloat> {
public:
  grayImage3d<GLfloat>() : grayBase3d<GLfloat>() {}
  grayImage3d<GLfloat>(ulong r, ulong c, ulong l) : grayBase3d<GLfloat>(r, c, l) {}
  void drawPixels(ulong l, bool inverseVideo = false) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_FLOAT, &d_voxel[l*d_rows*d_cols]);
  }
};

template<>
class grayImage3d<GLdouble> : public grayBase3d<GLdouble> {
public:
  grayImage3d<GLdouble>() : grayBase3d<GLdouble>() {}
  grayImage3d<GLdouble>(ulong r, ulong c, ulong l) : grayBase3d<GLdouble>(r, c, l) {}
  void drawPixels(ulong l, bool inverseVideo = false) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_DOUBLE, &d_voxel[l*d_rows*d_cols]);
  }
};

template<>
class grayImage3d<GLbyte> : public grayBase3d<GLbyte> {
public:
	//grayImage3d<GLbyte>() : grayBase3d<GLbyte>() {}
	//grayImage3d<GLbyte>(ulong r, ulong c, ulong l) : grayBase3d<GLbyte>(r, c, l) {}
  void drawPixels(ulong l, bool inverseVideo) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
	glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_BYTE, &d_voxel[l*d_rows*d_cols]);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
	glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
	}
};

template<>
class grayImage3d<GLuint> : public grayBase3d<GLuint> {
public:
	//grayImage3d<GLuint>() : grayBase3d<GLuint>() {}
	//grayImage3d<GLuint>(ulong r, ulong c, ulong l) : grayBase3d<GLuint>(r, c, l) {}
  void drawPixels(ulong l, bool inverseVideo) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_INT, &d_voxel[l*d_rows*d_cols]);
  }
};

template<>
class grayImage3d<GLushort> : public grayBase3d<GLushort> {
public:
	//grayImage3d<GLushort>() : grayBase3d<GLushort>() {}
	//grayImage3d<GLushort>(ulong r, ulong c, ulong l) : grayBase3d<GLushort>(r, c, l) {}
  void drawPixels(ulong l, bool inverseVideo) {
	assert(l < d_layers);
	setColorMap(inverseVideo);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
	glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_SHORT, &d_voxel[l*d_rows*d_cols]);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
	glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }
};

#endif /* __GRAYIMAGE3D__ */
