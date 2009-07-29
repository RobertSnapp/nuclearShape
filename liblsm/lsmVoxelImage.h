#ifndef __lsmVoxelImage_H__
#define __lsmVoxelImage_H__

//
//  lsmVoxelImage
//
//  Project liblsm
//
//  Created by Robert Snapp on 2007-04-28.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include "voxelImage.h"
#include "env.h"
#include "lsm.h"  // for type uint8, etc.
#include <vector>

class lsmVoxelImage : public voxelImage<uint8> {
private:
	LSM_File* d_lsm;
	std::vector<uint32> d_bandIndices;

public:
  lsmVoxelImage();

  // This constructor loads the indicated channel (or band) from
  // the lsm file indicated by the first argument.
  explicit lsmVoxelImage(char const *filename, int channel);
	
  //Thi constructor load all channels (or bands) from the
  // lsm file indicated by the first argument.
  explicit lsmVoxelImage(char const *filename);
	
  ~lsmVoxelImage();
	
  bool read(char const *filename);
  bool read(char const *filename, uint32 band);
  bool read(std::string const &filename);
  bool read(std::string const &filename, uint32 band);

  float getDx() { // voxel dx in microns.
    return static_cast<float>(d_lsm->getVoxelSizeDx());
  }
  float getDy() { // voxel dy in microns.
    return static_cast<float>(d_lsm->getVoxelSizeDy());
  }
  float getDz() { // voxel dz in microns.
    return static_cast<float>(d_lsm->getVoxelSizeDz());
  }

  void clear() {
    voxelImage<uint8>::clear();
    if (d_lsm) {
      delete d_lsm;
	  d_lsm = 0;
	}
	d_bandIndices.clear();
  }
	

};

#endif // __lsmVoxelImage_H__
