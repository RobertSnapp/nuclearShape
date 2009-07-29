#ifndef _IMAGE_COMPONENT_ARRAY_H_
#define _IMAGE_COMPONENT_ARRAY_H_
//
//  imageComponentArray.h
//  lsmModeler
//
//  Created by Robert Snapp on 2007-02-06.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include "env.h"  // for OpenGL and GLUT
#include "grayImage3d.h"
#include "lsmVoxelImage.h"
#include "nuclearSurface.h"
#include "convexhull.h"
#include "ntuple.h"
#include <utility>
#include <vector>

struct Subvolume;

enum HullType {VOXEL_CENTERS, VOXEL_CORNERS};


class ImageComponentArray {
private:
  grayImage3d<ulong> d_labels;   // changed to a three-dimensional image buffer.
  grayImage3d<GLubyte> d_colors;
  std::vector< std::pair<int,int> > d_labelCount;
  static const bool debug = true;

public:		
#ifdef COMMENT	
  ImageComponentArray(const lsmVoxelImage &image, 
							   ulong channel, 
							   Subvolume &subvol, 
							   GLubyte threshold = 0);
#endif

  ImageComponentArray() { };

//	explicit ImageComponentArray(const LSM_Image &image, GLubyte threshold = 0);


	ImageComponentArray(const ImageComponentArray &ica) {
		d_labels     = ica.d_labels;
		d_colors     = ica.d_colors;
		d_labelCount = ica.d_labelCount;
	}

	ImageComponentArray& operator=(ImageComponentArray &ica) {
		d_labels     = ica.d_labels;
		d_colors     = ica.d_colors;
		d_labelCount = ica.d_labelCount;
		return *this;
	}
  void computeComponents(const binaryImage3d &input);
  //  void computeComponents(const grayImage3d<GLubyte> &input);
  void computeConvexHull(ulong i, ConvexHull<int> &c, HullType mode=VOXEL_CENTERS);

  void clear() {
	d_labels.clear();
	d_colors.clear();
	d_labelCount.clear();
  }

  size_t size() {
	return d_labelCount.size();
  }
	
  void displayComponent(int i);
  void displayLayer(ulong l);
  void renderVoxel(int i, int j, int l);
  void floodSubstitute(int i, int j, int l, ulong oldValue, ulong newValue);

  void pushComponentVoxels(int rank, std::vector<ntuple<int, 3> > &points);
  void pushComponentShell(ulong rank, std::vector<ntuple<int, 3> > &points);

  ulong boundaryVoxelCount(int i) const;
  ulong clusterSize(ulong i) const  {
	assert(i < d_labelCount.size());
	return d_labelCount[i].second;}

  double volume(ulong i) const;

  void setDX(float dx) {
    d_labels.setDX(dx);
  }

  void setDY(float dy) {
    d_labels.setDY(dy);
  }

  void setDZ(float dz) {
    d_labels.setDZ(dz);
  }

  float getDX() {
    return d_labels.getDX();
  }

  float getDY() {
    return d_labels.getDY();
  }

  float getDZ() {
    return d_labels.getDZ();
  }

  ulong rows() { return d_labels.rows();}
  ulong cols() { return d_labels.cols();}
  ulong layers() { return d_labels.layers();}

  void makeColorLabels(size_t n);
};


#endif /* _IMAGE_COMPONENT_ARRAY_H_ */
