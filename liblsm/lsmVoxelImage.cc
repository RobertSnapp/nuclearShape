//
//  lsmVoxelImage.cc
//
//  Created by Robert Snapp on 2007-04-28.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include "lsm.h"
#include "lsmVoxelImage.h"

#include <algorithm>
#include <ios>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

lsmVoxelImage::lsmVoxelImage() : voxelImage<uint8>(), d_lsm(0) {}

lsmVoxelImage::lsmVoxelImage(char const *fileName, int bandIndex) : voxelImage<uint8>(), d_lsm(0) {
  bool status = this->read(fileName, bandIndex);
  if (! status) {
	throw std::ios_base::failure("lsmVoxelImage can't open file.");
  }
}

lsmVoxelImage::lsmVoxelImage(char const *fileName) : voxelImage<uint8>(), d_lsm(0) {
  bool status = this->read(fileName);
 if (! status) {
	throw std::ios_base::failure("lsmVoxelImage can't open file.");
  }
}

// lsmVoxelImage::lsmVoxelImage transfers the pixel data of the indicated single 
// channel image in an LSMfile, to the pixel buffer in a gray level image that 
// assumes one unsigned byte per pixel. 

bool lsmVoxelImage::read(std::string const &filename, uint32 bandIndex) {
  return this->read(filename.c_str(), bandIndex);
}


bool lsmVoxelImage::read(char const *fileName, uint32 bandIndex) {	
  bool status = true;
  uint32 l_rows;
  uint32 l_cols;
  uint32 l_bands;	
  uint32 layer = 0;
	
  try {
    d_lsm = (LSM_File *) new LSM_File(fileName);
  } 
  catch (std::ios_base::failure) {
    return false;
  }

  if (! d_lsm) {
    //	throw std::ios_base::failure("lsmVoxelImage:Can't open fileName.");
    // the indicated file is either not accessible, or is not an lsm file.
    return false;
  }
	
  LSM_File::iterator lsmEnd = d_lsm->end();
  for(LSM_File::iterator iter = d_lsm->begin(); iter != lsmEnd; ++iter) {
    if (! iter->isThumbnail() ) {
      if (++layer == 1) {
        l_rows  = iter->getLength();
        l_cols  = iter->getWidth();
        l_bands = iter->getSamplesPerPixel();
      } else {
        uint32 r = iter->getLength();
        uint32 c = iter->getWidth();
        uint32 b = iter->getSamplesPerPixel();
        if (r != l_rows || c != l_cols || b != l_bands) {
          std::cerr << "lsmVoxelImage: Warning: "
                    << " the dimensions of layer " << layer
                    << " differ from those in previous layer."
                    << std::endl;

          l_rows  = std::max(l_rows,  r);
          l_cols  = std::max(l_cols,  c);
          l_bands = std::max(l_bands, b);
        }
      }
    }
  }
		
  if (bandIndex < 0 || bandIndex >= l_bands) {
    throw std::out_of_range("lsmVoxelImage: bad bandIndex.");
  }
  d_bandIndices.push_back(bandIndex);
	
  voxelImage<uint8>::setRows(l_rows);
  voxelImage<uint8>::setCols(l_cols);
  voxelImage<uint8>::setLayers(layer);
  voxelImage<uint8>::setBands(1);
  voxelImage<uint8>::resize();
	
  double scale = 10.0/(l_rows*d_lsm->getVoxelSizeDx());
  voxelImage<uint8>::setDX(scale*d_lsm->getVoxelSizeDx());
  voxelImage<uint8>::setDY(scale*d_lsm->getVoxelSizeDy());
  voxelImage<uint8>::setDZ(scale*d_lsm->getVoxelSizeDz());
	
  uint32 l = 0;
  for(LSM_File::iterator iter = d_lsm->begin(); 
      iter != lsmEnd && l < layer; ++iter) {
    if (! iter->isThumbnail()) { // Skip over thumbnails.
      // Copy image data from the lsm file (*lsm) to the voxel image.
      char *c = (char *) voxelImage<uint8>::getVoxelOffset(0, 0, l++, 0);
      (void )iter->getPixels(bandIndex, c);
    }
  }

  return status;
}

// lsmVoxelImage::lsmVoxelImage transfers the pixel data of the indicated single 
// channel image in an LSMfile, to the pixel buffer in a gray level image that 
// assumes one unsigned byte per pixel. 

bool lsmVoxelImage::read(std::string const &fileName) {
  return  this->read(fileName.c_str());
}

bool lsmVoxelImage::read(char const *fileName) {	
  bool status = true;
  uint32 l_rows;
  uint32 l_cols;
  uint32 l_bands;	
  uint32 layer = 0;
	
//   d_lsm = (LSM_File *) new LSM_File(fileName);
//   if (! d_lsm) {
//     //	throw std::ios_base::failure("lsmVoxelImage:Can't open fileName.");
//     return false;
//   }

  try {
    d_lsm = (LSM_File *) new LSM_File(fileName);
  } 
  catch (std::ios_base::failure) {
    return false;
  }

  LSM_File::iterator lsmEnd = d_lsm->end();
  for(LSM_File::iterator iter = d_lsm->begin(); iter != lsmEnd; ++iter) {
    if (! iter->isThumbnail() ) {
      if (++layer == 1) {
        l_rows  = iter->getLength();
        l_cols  = iter->getWidth();
        l_bands = iter->getSamplesPerPixel();  // number of channels
      } else {
        uint32 r = iter->getLength();
        uint32 c = iter->getWidth();
        uint32 b = iter->getSamplesPerPixel();
        if (r != l_rows || c != l_cols || b != l_bands) {
          std::cerr << "lsmVoxelImage: Warning: "
                    << " the dimensions of layer " << layer
                    << " differ from previous layers."
                    << std::endl;
					
          l_rows =  std::max(l_rows,  r);
          l_cols =  std::max(l_cols,  c);
          l_bands = std::max(l_bands, b);
        }
      }
    }
  }
	
  for(uint32 bandIndex = 0; bandIndex < l_bands; bandIndex++)
    d_bandIndices.push_back(bandIndex);
		
  voxelImage<uint8>::setRows(l_rows);
  voxelImage<uint8>::setCols(l_cols);
  voxelImage<uint8>::setLayers(layer);
  voxelImage<uint8>::setBands(l_bands);
  voxelImage<uint8>::resize();
	
  voxelImage<uint8>::setDX(d_lsm->getVoxelSizeDx());
  voxelImage<uint8>::setDY(d_lsm->getVoxelSizeDy());
  voxelImage<uint8>::setDZ(d_lsm->getVoxelSizeDz());
	
  uint32 l = 0;
  for(LSM_File::iterator iter = d_lsm->begin(); 
      iter != lsmEnd && l < layer; ++iter) {
    if (! iter->isThumbnail()) { // skip thumbnails.
      for(uint16 b = 0; b < iter->getSamplesPerPixel(); b++) {
        // Copy image data from the lsm file (*lsm) to the voxel 
        // image.
        char *c = (char *) voxelImage<uint8>::getVoxelOffset(0, 0, l, b);
        (void) iter->getPixels(b, c);
      }
      l++;
    }
  }

  return status;
}

lsmVoxelImage::~lsmVoxelImage() {
  //  std::cerr << "deleting lsmVoxelImage." << std::endl;
  if (d_lsm) {
	delete d_lsm;
	d_lsm = 0;
  }
}
