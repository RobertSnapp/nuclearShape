//
//  imageComponentArray.cc
//  lsmModeler
//
//  Created by Robert Snapp on 2007-02-06.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include "nuclearSurface.h"
#include "imageComponentArray.h"
#include "grayImage3d.h"
#include "env.h"
#include "lsmVoxelImage.h"
#include "imageHistogram.h"
#include "ntuple.h"
#include <algorithm>
#include <stdexcept>

using namespace std;
const bool debug = true;

bool greaterFrequency(pair<int,int> lhs, pair<int,int> rhs) {
	return lhs.second > rhs.second;
}

// floodSubstitute uses a "flood fill" algorithm to efficiently substitute all
// occurences of the voxel label oldValue with newValue. If the voxel in row i,
// column j, layer l, is labeled with the oldValue, then the label is replaced
// with new label, and the algorithm is recursively applied to all neighboring
// voxels that occur before (i, j, l). Thus the algorithm is applied to the 26
// neighboring voxels. (Note, that if (i, j, l) falls on a boundary of the
// voxel image, then only a subset of neighbors are examined.)
void ImageComponentArray::floodSubstitute(int i, int j, int l, int oldValue, int newValue) {
	if (labels_.getVoxel(i, j, l) == oldValue) {
		labels_.setVoxel(i, j, l, newValue);
		
		for(int ll = max<int>(0, l-1); ll < min<int>(labels_.layers(), l+2); ll++)
			for(int ii = max<int>(0, i-1); ii < min<int>(labels_.rows(), i+2); ii++)
				for(int jj = max<int>(0, j-1); jj < min<int>(labels_.cols(), j+2); jj++) {
					floodSubstitute(ii, jj, ll, oldValue, newValue);
				}		
	}
	return;
}


// A component is defined as a group of adjacent voxels that all contain intensities that are
// greater than or equal to a given threshold. The following function, assigns a unique
// integer label to each component. After calling labelComponents, the output grayImage3d
// labels can be scanned to identify the number of components, and the size of each, that
// lie above a given threshold. If threshold = 0, then an adaptive threshold is computed.

// voxels below threshold are labeled with a 0; those above with a positive integer, such
// that each contiguous set of voxels (assuming a 26 point topology)

// TODO revise the adaptive threshold procedure.
// TODO paramaterize the topology to include 18 and 6 point topologies. 

ImageComponentArray::ImageComponentArray(	const lsmVoxelImage &input, 
											int channel,
                                            Subvolume &subvol,
											GLubyte threshold) {
	if (debug) {
		cerr << "ImageComponentArray(...):" << endl;
	}
	
	int nextLabel = 1;
	uint32 rows   = input.rows();
	uint32 cols   = input.cols();
	uint32 layers = input.layers();
	
	float dx = input.getDX();
	float dy = input.getDY();
	float dz = input.getDZ();
	float fthreshold = static_cast<float>(threshold);
	
	if (threshold == 0) {
		// select a threshold using the following adaptive algorithm
		ImageHistogram h(1);
		h.processImage(input, channel);
		fthreshold = static_cast<float>(h.getMean() + 1.5*h.getSDev());	
	}
	
	if(debug){
		cerr << "rows = " << rows << endl
		     << "cols = " << cols << endl
		     << "layers = " << layers << endl
		     << "dx = " << dx << " meters" << endl
			 << "dy = " << dy << " meters" << endl
			 << "dz = " << dz << " meters" << endl
	    	 <<"actual threshold = " << fthreshold << endl;
	}
	
	labels_.setDX(dx);
	labels_.setDY(dy);
	labels_.setDZ(dz);
	labels_.setRows(rows);
	labels_.setCols(cols);
	labels_.setLayers(layers);
	labels_.resize();
    
    int lmin = subvol.layer;
    int lmax = lmin + subvol.depth;
    int cmin = subvol.col;
    int cmax = cmin + subvol.width;
    int rmin = subvol.row;
    int rmax = rmin + subvol.height;

	for(int l = 0; l < subvol.depth; l++) {
		if (debug) cerr << "!"; 						// poor man's progress bar.
		for(int i = 0; i < subvol.height; i++) {
			if (debug && i % 128 == 0) cerr << ".";		// poor man's progress bar.
			for(int j = 0; j < subvol.width; j++) {
				if (input.getVoxelClip(i + rmin, j + cmin, l + lmin, channel) >= fthreshold) {
				//	Assign a label to the current voxel that matches the label of any contiguous
				//  neighbors, assuming a 26-point topology, of which at most 9 + 3 + 1 = 13
				//  adjacent neighbors have been examined.
				
					for(int ll = max(0, l-1); ll <= l; ll++) 
						for (int ii = max<int>(0, i-1); ii <  min<int>(0, i + l - ll + 1); ii++)
							for (int jj = max<int>(0, j-1); jj < min<int>(subvol.width, (ll < l || ii < i ? j + 2 : j)); jj++) {
								int neighbor = labels_.getVoxel(ii, jj, ll);
								if (neighbor > 0) {
									int center = labels_.getVoxel(i, j, l);
									if (center == 0) {
										labels_.setVoxel(i, j, l, static_cast<ulong>(neighbor));
									} else if (neighbor != center) {
										int oldValue = max<int>(neighbor, center);
										int newValue = min<int>(neighbor, center);
										
										floodSubstitute(i, j, l, oldValue, newValue);
									}
								}
							}					
					
					if (labels_.getVoxel(i, j, l) == 0) {                            
					 // no upper-left neighboring pixels are set. Select a new component label.
						labels_.setVoxel(i, j, l, static_cast<GLint>(nextLabel++));      
					}
				}
			}
		}
	}	
	
	if (debug) cerr << nextLabel << " labels, ";
	
	// Order the positively labeled components, according to their prominence. 
	labelCount_.resize(nextLabel);
	for(int k = 0; k < nextLabel; k++) {
		labelCount_[k].first = k;
	}
	
	for (int l = 0; l < layers; l++)
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++) {
				int k = labels_.getVoxel(i, j, l);
				if (k > 0) labelCount_[k].second += 1;
	}
	
	std::sort(labelCount_.begin(), labelCount_.end(), greaterFrequency);

	if (debug) {
		for(int i = 0; i < min(25, nextLabel); i++) {
			cerr << labelCount_[i].second << ", ";
		}
	}
	
	if (debug) cerr << "done!" << endl;
}

void ImageComponentArray::pushComponentVoxels(int rank, std::vector<ntuple<double, 3> > &points, 
	float dx, float dy, float dz) {
	if (rank < 0 || rank > labelCount_.size())
		throw std::out_of_range("pushComponentVoxels: bad value for rank.");
		
	int label = labelCount_[rank].first;
	for(int l = 0; l < labels_.layers(); ++l)
		for(int i = 0; i < labels_.rows(); ++i)
			for(int j = 0; j < labels_.cols(); ++j) {
				if (labels_.getVoxel(i, j, l) == label)
					points.push_back(ntuple<double, 3>(static_cast<double>(i*dx), 
												       static_cast<double>(j*dy), 
												       static_cast<double>(l*dz)));
			}
}

void ImageComponentArray::renderVoxel(int i, int j, int l) {
	float dx = labels_.getDX();
	float dy = labels_.getDY();
	float dz = labels_.getDZ();
	float x[2], y[2], z[2];
	
	for(int ii = 0; ii < 2; ii++) {
		x[ii] = (i + ii - 0.5)*dx;
		y[ii] = (j + ii - 0.5)*dy;
		z[ii] = (l + ii - 0.5)*dz;
	}
	
	glBegin(GL_QUADS);
		glNormal3f(0, 0, -1.);
		glVertex3f(x[0], y[0], z[0]);
		glVertex3f(x[0], y[1], z[0]);
		glVertex3f(x[1], y[1], z[0]);
		glVertex3f(x[1], y[0], z[0]);
		
		glNormal3f(0, -1., 0);
		glVertex3f(x[0], y[0], z[0]);
		glVertex3f(x[1], y[0], z[0]);
		glVertex3f(x[1], y[0], z[1]);
		glVertex3f(x[0], y[0], z[1]);
		
		glNormal3f(-1., 0., 0.);
		glVertex3f(x[0], y[0], z[0]);
		glVertex3f(x[0], y[0], z[1]);
		glVertex3f(x[0], y[1], z[1]);
		glVertex3f(x[0], y[1], z[0]);
		
		glNormal3f(0., 0., 1.);
		glVertex3f(x[0], y[0], z[1]);
		glVertex3f(x[1], y[0], z[1]);
		glVertex3f(x[1], y[1], z[1]);
		glVertex3f(x[0], y[1], z[1]);

		glNormal3f(0., 1., 0.);
		glVertex3f(x[0], y[1], z[0]);
		glVertex3f(x[0], y[1], z[1]);
		glVertex3f(x[1], y[1], z[1]);
		glVertex3f(x[1], y[1], z[0]);

		glNormal3f(1., 0., 0.);
		glVertex3f(x[1], y[0], z[0]);
		glVertex3f(x[1], y[1], z[0]);
		glVertex3f(x[1], y[1], z[1]);
		glVertex3f(x[1], y[0], z[1]);
	glEnd();						
}

void ImageComponentArray::displayComponent(int component) {
	if (debug) cerr << "displayComponent(" << component << ")... ";
	int targetLabel = labelCount_[component].first;
	int voxelCount = 0;
	for(int l = 0; l < labels_.layers(); l++)
		for(int i = 0; i < labels_.rows(); i++)
			for(int j = 0; j < labels_.cols(); j++) {
				if (targetLabel == labels_.getVoxel(i, j, l)) {
					renderVoxel(i, j, l);
					voxelCount++;
				}
			}
			
	if (debug) cerr << voxelCount << " voxels rendered ... done!" << endl;
	return;
}

ulong ImageComponentArray::boundaryVoxelCount(int component) const {
	ulong boundaryVoxels = 0;
	
	int layers = labels_.layers();
	int rows   = labels_.rows();
	int cols   = labels_.cols();
	
	int targetLabel = labelCount_[component].first;
	for(int l = 0; l < layers; l++)
		for(int i = 0; i < rows; i++)
			for(int j = 0; j < cols; j++) {
				if (labels_.getVoxelClip(i, j, l) == targetLabel) {
					bool isBoundaryVoxel = false;
					for(int ll = max<int>(0, l-1); ll < min<int>(layers, l+2); ll++)
						for(int ii = max<int>(0, i-1); ii < min<int>(rows, i+2); ii++)
							for(int jj = max<int>(0, j-1); jj < min<int>(cols, j+2); jj++) {
								if (labels_.getVoxelClip(ii, jj, ll) != targetLabel)
									isBoundaryVoxel = true;
							}
					if (isBoundaryVoxel) boundaryVoxels++;
				}
	}
	return boundaryVoxels;			
}
