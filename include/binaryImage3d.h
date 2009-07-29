/*****************************************************************************/
/* binaryImage3d.h 
/* 
/* Implements a three dimensional binary image.
/*
/* Robert R. Snapp  Copyright (C) 2007.
/*****************************************************************************/

#ifndef __BINARY_IMAGE_3D_H__
#define __BINARY_IMAGE_3D_H__

#include <algorithm> // min function
#include <cmath> // log function
#include <valarray>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cassert>
#include <typeinfo>
#include "voxelArray.h"
#include "env.h"

typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;


class binaryImage3d : public voxelArray {
protected:
  static const int bitsPerByte = 8;
  std::valarray<uchar> d_voxel;
  int d_bytesPerRow;
  int d_bytesPerLayer;
public:
	binaryImage3d(int rows, int cols, int layers) : voxelArray(rows, cols, layers, 1) {
		d_bytesPerRow = static_cast<int>(ceil(static_cast<double>(d_cols)/bitsPerByte));
		d_bytesPerLayer = d_rows*d_bytesPerRow;	
		d_voxel.resize(d_bytesPerLayer*d_layers);		
	}

 binaryImage3d() : voxelArray(), d_bytesPerRow(0), d_bytesPerLayer(0) {
		voxelArray::d_bands = 1;
	}
	
	~binaryImage3d() {}
		
	binaryImage3d(const binaryImage3d & p) : voxelArray(p) {
      d_bytesPerRow   = p.d_bytesPerRow;
      d_bytesPerLayer = p.d_bytesPerLayer;		
      d_voxel.resize(d_bytesPerLayer*d_layers);
      d_voxel = p.d_voxel;  // Deep copy.	
	}
		
	void clear() {
		voxelArray::clear();
		d_bytesPerRow   = 0;
		d_bytesPerLayer = 0;
		d_voxel.resize(0);
	}
	
	binaryImage3d& operator=(const binaryImage3d &p) {
		// guard against self-assignment
		if (this != &p) {
          d_bytesPerRow   = p.d_bytesPerRow;
          d_bytesPerLayer = p.d_bytesPerLayer;
          static_cast<voxelArray&>(*this) = static_cast<const voxelArray&>(p);
          d_voxel.resize(d_bytesPerLayer*d_layers);
          d_voxel = p.d_voxel;
		}
		return *this;
	}
	
  void setVoxel(ulong i, ulong j, ulong l, bool x) {
    assert(0 <= i); assert(i < d_rows);
    assert(0 <= j); assert(j < d_cols);
    assert(0 <= l); assert(l < d_layers);

    if (0 <= i && i < d_rows && 0 <= j && j < d_cols) {
      if (x)
        d_voxel[l*d_bytesPerLayer + i*d_bytesPerRow + j/bitsPerByte] |=  0x01 << (j % bitsPerByte);
      else
        d_voxel[l*d_bytesPerLayer + i*d_bytesPerRow + j/bitsPerByte] &= 0xff ^ (0x01 << (j % bitsPerByte));
    }
  }
    
  bool getVoxel(ulong i, ulong j, ulong l) const {
	assert(0 <= i); assert(i < d_rows);
    assert(0 <= j); assert(j < d_cols);
    assert(0 <= l); assert(l < d_layers);

    return d_voxel[l*d_bytesPerLayer + i*d_bytesPerRow + j/bitsPerByte] & (0x01 << (j % bitsPerByte));
  }

  bool getVoxelClip(ulong i, ulong j, ulong l) const {
    if (0 <= i && i < d_rows && 0 <= j && j < d_cols && 0 <= l && l < d_layers)
      return d_voxel[l*d_bytesPerLayer + i*d_bytesPerRow + j/bitsPerByte] & (0x01 << (j % bitsPerByte));
    else
      return 0;
  }
  
  void getVoxel(ulong i, ulong j, ulong l, bool &v) const {
    v = d_voxel[l*d_bytesPerLayer + i*d_bytesPerRow + j/bitsPerByte] & (0x01 << (j % bitsPerByte));
  }
    
 uchar* getVoxelOffset(ulong i, ulong j, ulong l) {
   return &d_voxel[l*d_rows*d_cols+ i*d_cols + j];
 }

  void resize() {
    d_bands = 1;
    d_bytesPerRow = static_cast<int>(ceil(static_cast<double>(d_cols)/bitsPerByte));
    d_bytesPerLayer = d_rows*d_bytesPerRow;
    d_voxel.resize(d_rows*d_cols*d_layers);
    d_bytesPerLayer = d_rows*d_bytesPerRow;
	}

  // drawPixels will display the image in the current OpenGL context, at the current raster position, using
  // the current zoom factor.
  void drawPixels(ulong l, bool inverseVideo = false) const {
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
		glDrawPixels(w, h, GL_COLOR_INDEX, GL_BITMAP, &d_voxel[l*d_bytesPerLayer]);
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
  }

};

 
#endif /* __BINARYIMAGE3D__ */
