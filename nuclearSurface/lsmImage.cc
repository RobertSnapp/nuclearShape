//
//  file: imageData.cc
//  project: lsmModeler
// 
//  manages the data acquired from an LSM image.
//
//  Created by Robert Snapp on 2007-01-12.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include "lsmImage.h"
#include "lsm.h"
#include "grayImage3d.h"
#include "lsmFileInfo.h"
#include "imageComponentArray.h"

using namespace std;

static const bool debug = true;

// loadImageData transfers the pixel data of the indicated single channel image in an LSM
// file, to the pixel buffer in a gray level image that assumes one unsigned byte per
// pixel. The number of bytes transferred is returned.

LSM_Image::LSM_Image(LSM_File* lsm, int channelIndex) : lsmFile_(lsm), channel_(channelIndex) {
	int nBytes = 0;
	
	info_ = (LSM_File_Info *) new LSM_File_Info(*lsm);
	info_->dump(std::cerr);
	int imageIndex = info_->findNextOriginalImage(0);
	int rows   = lsmFile_->getLength(imageIndex);
	std::cerr << "rows = " << rows << endl;
	int cols   = lsmFile_->getWidth(imageIndex);
	std::cerr << "cols = " << cols << endl;
	int layers = info_->numberOfOriginals();
	std::cerr << "layers = " << layers << endl;
	int depth  = lsmFile_->getSamplesPerPixel(imageIndex);
	std::cerr << "samples per pixel = " << depth << endl;
	
	if (channelIndex < 0 || channelIndex >= depth) {
		throw std::out_of_range("loadImage: bad channelIndex.");
	}
	
	image_.setRows(rows);
	image_.setCols(cols);
	image_.setLayers(layers);
	image_.resize();
	
	//double scale = 4.0/(rows*lsm->getVoxelSizeDx());
	image_.setDX(lsmFile_->getVoxelSizeDx());
	image_.setDY(lsmFile_->getVoxelSizeDy());
	image_.setDZ(lsmFile_->getVoxelSizeDz());
	
	for(int l = 0; l < layers; l++) {
		// load image data
		int newBytes = lsmFile_->getStripByteCount(imageIndex, channel_);
		if (newBytes != rows*cols) {
			std::cerr 	<< "Error: Data Inconsistency: in layer " << l 
						<< ", nBytes != rows*cols." << endl;
			throw bad_format();
		}
	
		// Copy image data from the lsm file (*lsm) to the grayLevel image (inputImage)
		char *c = (char *) image_.getVoxelOffset(0, 0, l);
		(void) lsm->getPixels(imageIndex, channel_, c);
		nBytes += newBytes;
		imageIndex = info_->findNextOriginalImage(imageIndex + 1);
	}
	
	std::cerr << "LSM_Image contstuctor data field test:" << endl;
	std::cerr << "image_.rows() = " << image_.rows() << endl;
	std::cerr << "image_.cols() = " << image_.cols() << endl;
	std::cerr << "image_.layers() = " << image_.layers() << endl;
	std::cerr << "bytes loaded = " << nBytes << endl << endl;
}

#ifdef COMMENT
LSM_Image::LSM_Image(const char *fileName, int channelIndex) : channel_(channelIndex)  {	
	lsmFile_ = (LSM_File *) new LSM_File(fileName);
	if (! lsmFile_) {
		cerr << "ERROR: Could not open file " << fileName << endl;
		throw bad_filename();
	}

	LSM_Image(lsmFile_, channel_);	
}
#endif

#ifdef COMMENT
void LSM_Image::applyThreshold(int threshold) {
	if (debug) {
		cerr << "Plotting lsm data with threshold = " << threshold << " ... ";
	}
	
	// int hits = drawBoundaryBricks(imageData, threshold);
	components_ = ImageComponentArray(image_, static_cast<GLubyte>(threshold));
	
	// display the prominent component
	components_.displayComponent(0);
	
	ulong bvc = components_.boundaryVoxelCount(0);
	ulong cs  = components_.clusterSize(0);
	cout << "Boundary Voxels = " << bvc << endl;
	cout << "Cluster Size    = " << cs << endl;
	cout << "Ratio           = " << static_cast<double>(bvc)/static_cast<double>(cs) << endl;
	
	// debug(cerr << hits << " hits." << endl);
	return;
}
#endif
