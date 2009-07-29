/*****************************************************************************
 * grayImage.h 
 * 
 * Robert R. Snapp  Copyright (C) 2006.
 *
 * Defines a class for storing and accessing a monochrome (gray level)
 * raster image. C++ templates are used to facilitate images that use
 * a variety of pixel types, e.g. bytes, short ints, or long ints, for
 * different pixel depths.
 *****************************************************************************/

#ifndef __GRAYIMAGE_H__
#define __GRAYIMAGE_H__

#include <algorithm> // min function
#include <cmath> // log function
#include <valarray>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cassert>
#include <typeinfo>
#include "pixelArray.h"
#include "env.h"


typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;

template<typename T> 
class grayBase : public pixelArray {
protected:
	std::valarray<T> d_pixel;
public:
	grayBase(ulong r, ulong c) : pixelArray(r, c, 1) {
		assert(r > 0);
		assert(c > 0);
		pixelArray::d_bytesPerRow = c * sizeof(T);		
		d_pixel.resize(r * c);		
	}

	grayBase() : pixelArray() {
		pixelArray::d_bands = 1;
		d_bytesPerRow = 0;
		d_pixel = 0;
	}
	
	virtual ~grayBase() {}
		
	grayBase(const grayBase & p) {		
		d_pixel.resize(d_rows*d_cols);
		d_pixel = p.d_pixel;  // Deep copy.
		d_bytesPerRow = d_cols*sizeof(T);
	}
		
	void clear() {
		pixelArray::clear();
		d_pixel.resize(0);
	}
	
	grayBase& operator=(const grayBase &p) {
      if (this != &p) {
		d_pixel.resize(d_rows*d_cols);
		d_pixel = p.d_pixel;
      }
      return *this;
	}
	
	T& graylevel(ulong i, ulong j) {
		assert(0 <= i); assert(i < d_rows);
		assert(0 <= j); assert(j < d_cols);
		return d_pixel[(i*d_cols + j)];
	}
			
    void setPixel(ulong i, ulong j, T v) {
    	assert(0 <= i); assert(i < d_rows);
		assert(0 <= j); assert(j < d_cols);
    	d_pixel[i*d_cols + j] = v;
    }
    
    void getPixel(ulong i, ulong j, T &v) {
    	v = d_pixel[i*d_cols + j];
    }
    
    T getPixel(ulong i, ulong j) {
    	if (0 <= i && i < d_rows && 0 <= j && j < d_cols) {
    		return d_pixel[i*d_cols + j];
    	} else {
    		return 0;
    	}
    }
    
 	T* getPixelOffset(ulong i, ulong j) {
		return &d_pixel[i*d_cols + j];
	}

    void setRows(ulong r) {
    	assert(r > 0); 
    	d_rows = r;
    }
    void setCols(ulong c) {
    	assert(c > 0); 
    	d_cols = c;
    }
   
    void resize() {
    	d_bands = 1;
    	d_pixel.resize(d_rows*d_cols);
		d_bytesPerRow = d_cols*sizeof(T);
	}
    
    double intensity() {
    	double val=0;
    	for(ulong i=0; i < d_rows; i++)
    		for(ulong j=0; j < d_cols; j++) {
    			val += static_cast<double>(d_pixel[i*d_cols + j]);
    	}
    	return val;
    }
    
    double entropy(ulong dr, ulong dc) {
    	assert(dr > 0);
    	assert(dc > 0);
    	double h=0;
    	unsigned long total = 0;
    	
    	for (ulong i = 0; i < d_rows; i++)
    		for (ulong j = 0; j < d_cols; j++) {
    			ulong cellIntensity = 0;
    			for (ulong k = i-dr/2; k < i+dr/2; k++)
    				for (ulong l = j-dc/2; l < j+dc/2; l++) {
    					cellIntensity += getPixel(k, l);
    			}
    			if (cellIntensity > 0) {
    				h -=  cellIntensity * std::log(static_cast<double>(cellIntensity));
    				total += cellIntensity;
    			}
    	}
    	return std::log(static_cast<double>(total)) + h/total;
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
					   		
  // virtual void drawPixels() = 0; 
};


template<typename T>
class grayImage : public grayBase<T> {
public:
	grayImage<T>() : grayBase<T>() {}
	grayImage<T>(ulong r, ulong c) : grayBase<T>(r, c) {}
};

template<>
class grayImage<GLubyte> : public grayBase<GLubyte> {
public:
  grayImage<GLubyte>() : grayBase<GLubyte>() {}
  grayImage<GLubyte>(const grayImage<GLubyte> &i) : grayBase<GLubyte>(static_cast<grayBase<GLubyte> >(i)) {}
  grayImage<GLubyte>(ulong r, ulong c) : grayBase<GLubyte>(r, c) {}

  grayImage<GLubyte>& operator=(const grayImage<GLubyte> &i) {
    if (this != &i) {
      grayBase<GLubyte>::operator=(static_cast<grayBase<GLubyte> >(i));
    }
    return *this;
  }

  void drawPixels(bool inverseVideo = false) {
	setColorMap(inverseVideo);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glDrawPixels(static_cast<GLsizei>(d_cols), 
                 static_cast<GLsizei>(d_rows), GL_LUMINANCE, GL_UNSIGNED_BYTE, &d_pixel[0]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }

  void rgbPixels() {
	setColorMap(false);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glDrawPixels(static_cast<GLsizei>(d_cols), 
                 static_cast<GLsizei>(d_rows), GL_RGB, GL_UNSIGNED_BYTE_3_3_2, &d_pixel[0]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }
};

template<>
class grayImage<GLshort> : public grayBase<GLshort> {
public:
  grayImage<GLshort>() : grayBase<GLshort>() {}
  grayImage<GLshort>(const grayImage<GLshort> &i) : grayBase<GLshort>(static_cast<grayBase<GLshort> >(i)) {}
  grayImage<GLshort>(ulong r, ulong c) : grayBase<GLshort>(r, c) {}

  grayImage<GLshort>& operator=(const grayImage<GLshort> &i) {
    if (this != &i) {
      grayBase<GLshort>::operator=(static_cast<grayBase<GLshort> >(i));
    }
    return *this;
  }
  void drawPixels(bool inverseVideo = false) {
	setColorMap(inverseVideo);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_SHORT, &d_pixel[0]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  }
};

template<>
class grayImage<GLint> : public grayBase<GLint> {
public:
  grayImage<GLint>() : grayBase<GLint>() {}
  grayImage<GLint>(ulong r, ulong c) : grayBase<GLint>(r, c) {}
  void drawPixels(bool inverseVideo = false) {
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_INT, &d_pixel[0]);
  }

  void mapPixels() {
    glDrawPixels(static_cast<GLsizei>(d_cols), 
                 static_cast<GLsizei>(d_rows), GL_COLOR_INDEX, GL_INT, &d_pixel[0]);
  }
};

template<>
class grayImage<GLfloat> : public grayBase<GLfloat> {
public:
  grayImage<GLfloat>() : grayBase<GLfloat>() {}
  grayImage<GLfloat>(ulong r, ulong c) : grayBase<GLfloat>(r, c) {}
  void drawPixels(bool inverseVideo = false) {
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_FLOAT, &d_pixel[0]);
  }
};

template<>
class grayImage<GLdouble> : public grayBase<GLdouble> {
public:
  grayImage<GLdouble>() : grayBase<GLdouble>() {}
  grayImage<GLdouble>(ulong r, ulong c) : grayBase<GLdouble>(r, c) {}
  void drawPixels(bool inverseVideo = false) {
	setColorMap(inverseVideo);
	glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_DOUBLE, &d_pixel[0]);
  }
};

template<>
class grayImage<GLbyte> : public grayBase<GLbyte> {
public:
	//grayImage<GLbyte>() : grayBase<GLbyte>() {}
	//grayImage<GLbyte>(ulong r, ulong c) : grayBase<GLbyte>(r, c) {}
	void drawPixels(bool inverseVideo) {
	  setColorMap(inverseVideo);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
      glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
      glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_BYTE, &d_pixel[0]);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
      glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
	}
};

template<>
class grayImage<GLuint> : public grayBase<GLuint> {
public:
	//grayImage<GLuint>() : grayBase<GLuint>() {}
	//grayImage<GLuint>(ulong r, ulong c) : grayBase<GLuint>(r, c) {}
	void drawPixels(bool inverseVideo) {
	  	setColorMap(inverseVideo);
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_INT, &d_pixel[0]);
	}
};

template<>
class grayImage<GLushort> : public grayBase<GLushort> {
public:
	//grayImage<GLushort>() : grayBase<GLushort>() {}
	//grayImage<GLushort>(ulong r, ulong c) : grayBase<GLushort>(r, c) {}
	void drawPixels(bool inverseVideo) {
	  setColorMap(inverseVideo);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 
      glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
      glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_SHORT, &d_pixel[0]);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 4); // restore default value
      glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
	}
};



#endif /* __GRAYIMAGE_H__ */
