//
//  morphology.cc
//  nuclearSlice
//
// Created by Robert R. Snapp on 2007-07-24
// Copyright © Robert R. Snapp 2007. All rights reserved.
//
// This file supports fundamental morphological operations for image analysis.
//


#include "binaryImage.h"
#include "project.h"

using namespace std;

binaryImage		   roundedSquare(4, 4);
binaryImage		   roundedSquareMask(4, 4);
binaryImage		   cross(3,3);
binaryImage        mediumSquare(3, 3);
binaryImage		   mediumSquareMask(3, 3);
binaryImage        smallSquare(2, 2);
binaryImage		   smallSquareMask(2, 2);


void initKernels() {
  if (verbosity > 10) {
	cerr << "Calling initKernels()..." << endl;
  }

  for(ulong i = 0; i < mediumSquare.rows(); i++) {
	for (ulong j = 0; j < mediumSquare.cols(); j++) {
	  mediumSquare.setPixel(i, j, 1);
	  mediumSquareMask.setPixel(i, j, 1);
	}
  }
	
  for(ulong i = 0; i < smallSquare.rows(); i++) {
	for (ulong j = 0; j < smallSquare.cols(); j++) {
	  smallSquare.setPixel(i, j, 1);
	  smallSquareMask.setPixel(i, j, 1);
	}
  }
	
  for(ulong i = 0; i < cross.rows(); i++) {
	for (ulong j = 0; j < cross.cols(); j++) {
	  bool value = (i == 1 || j == 1);
	  cross.setPixel(i, j, value);
	}
  }
		
  for(ulong i = 0; i < roundedSquare.rows(); i++) {
	for (ulong j = 0; j < roundedSquare.cols(); j++) {
	  roundedSquare.setPixel(i, j, 1);
	  bool maskVal = (i == 1 || i == 2 || j == 1 || j == 2);
	  roundedSquareMask.setPixel(i, j, maskVal);
	}
  }
  if (verbosity > 10) {
	cerr << "Returning from initKernels()." << endl;
  }
  return;
}


/* applyMorphAnd applies a morphological operation using the given mask and kernel. The suffix "And"
 * implies that the pixel at (i, j) is set if and only if every bit in the kernel that falls in the
 * support of the mask, matches every corresponding bit in the neigborhood of (i, j) in the image. */
void applyMorphAnd(binaryImage &inImg, binaryImage &morphkernel, binaryImage &mask, binaryImage &outImg) {
  if (verbosity > 10) {
	cerr << "Calling applyMorphAnd(binaryImage, binaryImage, binaryImage, binaryImage)..." << endl;
  }
  outImg.setRows(inImg.rows());
  outImg.setCols(inImg.cols());
  outImg.resize();
	
  ulong mrows = morphkernel.rows();
  ulong mcols = morphkernel.cols();
  
  for (ulong i = 0; i < outImg.rows(); i++) {
	for (ulong j = 0; j < outImg.cols(); j++) {
	  bool val = 1;
	  for (ulong k = 0; k < mrows; k++) {
		for (ulong l = 0; l < mcols; l++) {
		  if (mask.getPixel(k, l)) {
			val &= inImg.getPixel(i+k-mrows/2, j+l-mcols/2) == morphkernel.getPixel(k, l);
		  }}}
	  outImg.setPixel(i, j, val);	
	}
  }
   if (verbosity > 10) {
	cerr << "Returning from applyMorphAnd(binaryImage, binaryImage, binaryImage, binaryImage)." << endl;
  }
  return;
}


void applyMorphAnd(grayImage<GLubyte> &inImg, binaryImage &morphkernel, binaryImage &mask, grayImage<GLubyte> &outImg) {
  if (verbosity > 10) {
	cerr << "Calling applyMorphAnd(grayImage, grayImage, grayImage, grayImage)..." << endl;
  }
  outImg.setRows(inImg.rows());
  outImg.setCols(inImg.cols());
  outImg.resize();
	
  ulong mrows = morphkernel.rows();
  ulong mcols = morphkernel.cols();
	
  for (ulong i = 0; i < outImg.rows(); i++) {
	for (ulong j = 0; j < outImg.cols(); j++) {
	  bool val = 1;
	  for (ulong k = 0; k < mrows; k++) {
		for (ulong l = 0; l < mcols; l++) {
		  if (mask.getPixel(k, l)) {
			val &= inImg.getPixel(i+k-mrows/2, j+l-mcols/2) & 1 == morphkernel.getPixel(k, l);
		  }}}
	  if (val > 0) {
		outImg.setPixel(i, j, 0xff);
	  } else {
		outImg.setPixel(i, j, 0x00);
	  }
	}
  }
  if (verbosity > 10) {
	cerr << "Returning from applyMorphAnd(grayImage, grayImage, grayImage, grayImage)." << endl;
  }
  return;
}

/* applyMorphOr applies a morphologica operation using the given mask and kernel. The suffix "Or"
 * indicates that the pixel at (i, j) is set if any bit in the kernel, within the support of the mask,
 * matches a corresponding bit in the neigbhorhood of (i, j) in the image. */
void applyMorphOr(binaryImage &inImg, binaryImage &morphkernel, binaryImage &mask, binaryImage &outImg) {
  if (verbosity > 10) {
	cerr << "Calling applyMorpOr(binaryImage, binaryImage, binaryImage, binaryImage)..." << endl;
  }
  outImg.setRows(inImg.rows());
  outImg.setCols(inImg.cols());
  outImg.resize();
	
  ulong mrows = morphkernel.rows();
  ulong mcols = morphkernel.cols();
	
  for (ulong i = 0; i < outImg.rows(); i++) {
	for (ulong j = 0; j < outImg.cols(); j++) {
	  bool val = 0;
	  for (ulong k = 0; k < mrows; k++) {
		for (ulong l = 0; l < mcols; l++) {
		  if (mask.getPixel(k, l)) {
			val |= inImg.getPixel(i+k-mrows/2, j+l-mcols/2) == morphkernel.getPixel(k, l);
		  }}}
	  outImg.setPixel(i, j, val);	
	}
  }

  if (verbosity > 10) {
	cerr << "Returning from applyMorphOr()." << endl;
  }
  return;
}

void applyMorphOr(grayImage<GLubyte> &inImg, binaryImage &morphkernel, binaryImage &mask, grayImage<GLubyte> &outImg) {
  if (verbosity > 10) {
	cerr << "Calling applyMorphOr(grayImage, grayImage, grayImage, grayImage)..." << endl;
  }
  outImg.setRows(inImg.rows());
  outImg.setCols(inImg.cols());
  outImg.resize();
	
  ulong mrows = morphkernel.rows();
  ulong mcols = morphkernel.cols();
	
  for (ulong i = 0; i < outImg.rows(); i++) {
	for (ulong j = 0; j < outImg.cols(); j++) {
	  bool val = 0;
	  for (ulong k = 0; k < mrows; k++) {
		for (ulong l = 0; l < mcols; l++) {
		  if (mask.getPixel(k, l)) {
			val |= inImg.getPixel(i+k-mrows/2, j+l-mcols/2) & 1 == morphkernel.getPixel(k, l);
		  }}}
	  if (val > 0) {
		outImg.setPixel(i, j, 0xff);
	  } else {
		outImg.setPixel(i, j, 0x00);
	  }
	}
  }
  if (verbosity >10) {
	cerr << "Returning from applyMorphOr()." << endl;
  }
  return;
}


void erosionOp(binaryImage &inImg, binaryImage &outImg) {
  applyMorphAnd(inImg, smallSquare, smallSquareMask, outImg);
}

void erosionOp(grayImage<GLubyte> &inImg, grayImage<GLubyte> &outImg) {
  applyMorphAnd(inImg, smallSquare, smallSquareMask, outImg);
}


void dialationOp(binaryImage &inImg, binaryImage &outImg) {
  applyMorphOr(inImg, mediumSquare, mediumSquareMask, outImg);
}

void dialationOp(grayImage<GLubyte> &inImg, grayImage<GLubyte> &outImg) {
  applyMorphOr(inImg, mediumSquare, mediumSquareMask, outImg);
}

void closureOp() {
  applyMorphOr(afterThreshold, mediumSquare, mediumSquareMask, processedImage[0]);
  applyMorphAnd(processedImage[0], mediumSquare, mediumSquareMask, processedImage[1]);
}
