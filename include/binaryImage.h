/*****************************************************************************
 * binaryImage.h 
 * 
 * Robert R. Snapp Copyright (C) 2006.
 *****************************************************************************/

#ifndef __BINARY_IMAGE_H__
#define __BINARY_IMAGE_H__

#include <valarray>
#include <iostream>
#include <ostream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <typeinfo>
#include "pixelArray.h"
#include "env.h"

typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;

class binaryImage : public pixelArray {
private:
  static const int d_bitsPerByte = 8;
  GLubyte *d_pixel;

public:
  // Create binaryImage object with r rows and c columns
  binaryImage(int r, int c) : pixelArray(r, c, 1) {
    assert(r > 0);
    assert(c > 0);
    pixelArray::d_bytesPerRow = c / d_bitsPerByte + (c % d_bitsPerByte != 0 ? 1: 0);
	  // static_cast<int>(ceil(static_cast<double>(d_cols)/d_bitsPerByte));		
    // d_pixel.resize(d_rows*pixelArray::d_bytesPerRow);
    d_pixel = new GLubyte[r * pixelArray::d_bytesPerRow];
  }

 binaryImage() : pixelArray(0, 0, 1), d_pixel(0) {}
	
 ~binaryImage() {
    delete [] d_pixel;
  }
		
  binaryImage(const binaryImage & p) : pixelArray(p) {		
    //    d_pixel.resize(d_rows*d_bytesPerRow);
    // d_pixel = p.d_pixel;  // Deep copy.

    delete [] d_pixel;
    int pixelSize = d_rows*pixelArray::d_bytesPerRow;
    d_pixel = new GLubyte[pixelSize];
    for(int i = 0; i < pixelSize; i++) {
      d_pixel[i] = p.d_pixel[i];
    }
  }
		
  binaryImage& operator=(const binaryImage &p) {
    // guard against self-assignment
    if (this != &p) {
      this->pixelArray::operator=(p);
      delete [] d_pixel;
      int pixelSize = d_rows*pixelArray::d_bytesPerRow;
      d_pixel = new GLubyte[pixelSize];
      for(int i = 0; i < pixelSize; i++) {
        d_pixel[i] = p.d_pixel[i];
      }
      // d_pixel.resize(d_rows*d_bytesPerRow);
      // d_pixel = p.d_pixel;
    }
    return *this;
  }
	
  void clear() {
    pixelArray::clear();
    //    d_pixel.resize(0);
    delete [] d_pixel;
    d_pixel = 0;
  }
			
  void setPixel(ulong i, ulong j, bool x) {
    if (0 <= i && i < d_rows && 0 <= j && j < d_cols) {
      if (x)
        d_pixel[i*d_bytesPerRow + j/d_bitsPerByte] |=  0x01 << (j % d_bitsPerByte);
      else
        d_pixel[i*d_bytesPerRow + j/d_bitsPerByte] &= 0xff ^ (0x01 << (j % d_bitsPerByte));
    }
  }
    
  bool getPixel(ulong i, ulong j) const {
    if (0 <= i && i < d_rows && 0 <= j && j < d_cols)
      return d_pixel[i*d_bytesPerRow + j/d_bitsPerByte] & (0x01 << (j % d_bitsPerByte));
    else
      return 0;
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
    d_bytesPerRow = d_cols / d_bitsPerByte + (d_cols % d_bitsPerByte != 0 ? 1 : 0);
	  // static_cast<int>( ceil( static_cast<double>(d_cols)/d_bitsPerByte ) );
    //    d_pixel.resize(d_rows*d_bytesPerRow);
    delete [] d_pixel;
    d_pixel = new GLubyte[d_rows*d_bytesPerRow];
  }
  
  // drawPixels will display the image in the current OpenGL context, at the current raster position, using
  // the current zoom factor.
  void drawPixels(bool inverseVideo = false) const {
	GLfloat pos[4];  // Used to draw the image in discrete tiles.

	GLfloat map[2] = {0.0, 1.0};
	GLfloat zoomX;
	GLfloat zoomY;
	GLfloat x0;
	GLfloat y0;

	if (inverseVideo) {
	  map[0] = 1.0;
	  map[1] = 0.0;
	}

	glGetFloatv(GL_CURRENT_RASTER_POSITION, pos);
	glGetFloatv(GL_ZOOM_X, &zoomX);
	glGetFloatv(GL_ZOOM_Y, &zoomY);

	x0 = pos[0]/std::abs(zoomX);  // Initial raster coordinates
	y0 = pos[1]/std::abs(zoomY);  

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_UNPACK_LSB_FIRST, GL_TRUE);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_R, 2, map);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_G, 2, map);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_B, 2, map);

#ifdef NAIVE_IMPLEMENTAION
	// This naive implementation uses a single call to glDrawPixels. Small images (e.g., 256 x 256) appear 
	// correctly on my MacBook Pro. However, larger images (e.g., 1024 x 1024) are corrupted by a strange
	// artifact: the right half of the binary image is replaced with a copy of the left half, i.e., the first
	// 512 columns are wrapped. 
    glDrawPixels(static_cast<GLsizei>(d_cols), static_cast<GLsizei>(d_rows), 
                 GL_COLOR_INDEX, GL_BITMAP, d_pixel); // &d_pixel[0]);
#endif

	// Draw the pixels as a rectangular array of tiles. If the dimensions of each tile is sufficiently small
	// then the artifact noted above is avoided.

	GLsizei const tile_dc = 512;  // tile width
	GLsizei const tile_dr = 512;  // tile height
	GLsizei tile_cols = d_cols/tile_dc + (d_cols % tile_dc > 0 ? 1 : 0); // number of tile cols (discrete ceiling)
	GLsizei tile_rows = d_rows/tile_dr + (d_rows % tile_dr > 0 ? 1 : 0); // number of tile rows (discrete ceiling)
	GLsizei margin_dc = d_cols - (tile_cols - 1)*tile_dc;  // The number of columns in the rightmost
	GLsizei margin_dr = d_rows - (tile_rows - 1)*tile_dr;
	GLfloat dx = tile_dc * (zoomX > 0 ? 1.0 : -1.0);  // Raster coordinate horizontal offset per tile
	GLfloat dy = tile_dr * (zoomY > 0 ? 1.0 : -1.0);  // Raster coordinate vertical offset per tile
	glPixelStorei(GL_UNPACK_ROW_LENGTH, d_cols);
	for (GLsizei i = 0; i < tile_rows; i++) {
	  glPixelStorei(GL_UNPACK_SKIP_ROWS, i*tile_dr);
	  GLsizei h = (i < tile_rows - 1 ? tile_dr : margin_dr);
	  for (GLsizei j = 0; j < tile_cols; j++) {
		glPixelStorei(GL_UNPACK_SKIP_PIXELS, j*tile_dc);
		GLsizei w = (j < tile_cols - 1 ? tile_dc : margin_dc);
		glRasterPos2f(x0 + j*dx, y0 + i*dy);         // World coordinates of the tile's initial corner.
		glDrawPixels(w, h, GL_COLOR_INDEX, GL_BITMAP, &d_pixel[0]);
	  }
	}
		
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
    glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
    glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);

	GLenum errCode;
	const GLubyte *errString;
	while ((errCode = glGetError()) != GL_NO_ERROR) {
		errString = gluErrorString(errCode);
		std::cerr << "OpenGL Error detected in binaryImage::drawPixels(): " << errString << std::endl;
	}

  //  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);  // restore default value
  //  glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);  // restore default value
  }
    
  int count(bool val) const {
    int k = 0;
    for (ulong i = 0; i < d_rows; i++) {
      for (ulong j = 0; j < d_cols; j++) {
        if (getPixel(i, j) == val) k++;
      }
    }
    return k;
  }
    
  double energyRatio() const {
    long int activeBoundaries = 0;
    long int activePixels = 0;
    	
    for (ulong i = 0; i < d_rows; i++)
      for (ulong j = 0; j < d_cols; j++) {
        if (getPixel(i, j)) {
          activePixels++;
          if (!getPixel(i-1, j-1)) activeBoundaries++;
          if (!getPixel(i-1, j  )) activeBoundaries++;
          if (!getPixel(i-1, j+1)) activeBoundaries++;
    				
          if (!getPixel(i,   j-1)) activeBoundaries++;
          if (!getPixel(i,   j+1)) activeBoundaries++;
    				
          if (!getPixel(i+1, j-1)) activeBoundaries++;
          if (!getPixel(i+1, j  )) activeBoundaries++;
          if (!getPixel(i+1, j+1)) activeBoundaries++;
        }
      }
    return static_cast<double>(activeBoundaries)/activePixels;
  }

  void describe(std::ostream &os) const {
    this->pixelArray::describe(os);
    os << "d_bytesPerRow = " << d_bytesPerRow << ", "
       << "d_bitsPerByte = " << d_bitsPerByte << ", "
       << "d_pixel.size() = " << d_rows*d_bytesPerRow
       << std::endl;
  }
};



#endif // __BINARY_IMAGE_H__
