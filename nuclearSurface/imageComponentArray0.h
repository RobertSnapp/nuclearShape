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
#include "ntuple.h"
#include <utility>
#include <vector>

typedef unsigned long int ulong;

class ImageComponentArray {
private:
	grayImage3d<ulong> labels_;
	std::vector< std::pair<int,int> > labelCount_;
	
public:		
	explicit ImageComponentArray(const lsmVoxelImage &image, int channel, GLubyte threshold = 0);
//	explicit ImageComponentArray(const LSM_Image &image, GLubyte threshold = 0);
	
	ImageComponentArray(const ImageComponentArray &ica) {
		labels_     = ica.labels_;
		labelCount_ = ica.labelCount_;
	}

	ImageComponentArray& operator=(ImageComponentArray &ica) {
		labels_     = ica.labels_;
		labelCount_ = ica.labelCount_;
		return *this;
	}
	
	void displayComponent(int i);
	void renderVoxel(int i, int j, int l);
	void floodSubstitute(int i, int j, int l, int oldValue, int newValue);
	void pushComponentVoxels(int rank, std::vector<ntuple<double, 3> > &points,
		float dx, float dy, float dz);
	ulong boundaryVoxelCount(int i) const;
	ulong clusterSize(int i) const  {return labelCount_[i].second;}
};




#endif /* _IMAGE_COMPONENT_ARRAY_H_ */
