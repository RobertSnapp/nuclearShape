//
//  morphology.cc
//  nuclearSurface
//
// Created by Robert R. Snapp on 2007-07-24
// Copyright © Robert R. Snapp 2007. All rights reserved.
//
// This file supports fundamental three-dimensional morphological operations for image analysis.
//

#include <algorithm>
#include "binaryImage3d.h"
#include "grayImage3d.h"
#include "nuclearSurface.h"

using namespace std;

binaryImage3d roundedCube(4, 4, 4);
binaryImage3d roundedCubeMask(4, 4, 4);
binaryImage3d cross(3,3, 3);
binaryImage3d mediumCube(3, 3, 3);
binaryImage3d mediumCubeMask(3, 3, 3);
binaryImage3d smallCube(2, 2, 2);
binaryImage3d smallCubeMask(2, 2, 2);


void initKernels() {
  for(ulong l = 0; l < mediumCube.layers(); l++) {
	for(ulong i = 0; i < mediumCube.rows(); i++) {
      for (ulong j = 0; j < mediumCube.cols(); j++) { 
        mediumCube.setVoxel(i, j, l, 1);
        mediumCubeMask.setVoxel(i, j, l, 1);
      }
	}
  }
	
  for(ulong l = 0; l < smallCube.layers(); l++) {
	for(ulong i = 0; i < smallCube.rows(); i++) {
      for (ulong j = 0; j < smallCube.cols(); j++) {
        smallCube.setVoxel(i, j, l, 1);
        smallCubeMask.setVoxel(i, j, l, 1);
      }
	}
  }
	

  // Voxels in the three dimensional cross have at least two coordinates equal to 1. 
  for(ulong l = 0; l < cross.layers(); l++) {
	for(ulong i = 0; i < cross.rows(); i++) {
      for (ulong j = 0; j < cross.cols(); j++) {
        bool value = (i == 1 && j == 1) || (i ==1 && l == 1) || (j ==1 && l == 1);
        cross.setVoxel(i, j, l, value);
      }
	}
  }

  // Voxels in the rounded cube have at least one coordinate equal to 1.
  for(ulong l = 0; l < roundedCube.layers(); l++) {
	for(ulong i = 0; i < roundedCube.rows(); i++) {
      for (ulong j = 0; j < roundedCube.cols(); j++) {	
        bool maskVal = (i == 1) || (i == 2) || (j == 1) || (j ==2) || (l == 1) || (l == 2) ;
        roundedCube.setVoxel(i, j, l, 1);             
        roundedCubeMask.setVoxel(i, j, l, maskVal);
      }
	}
  }
	return;
}


/* applyMorphAnd applies a morphological operation using the given mask and kernel. The suffix "And"
 * implies that the pixel at (i, j) is set if and only if every bit in the kernel that falls in the
 * support of the mask, matches every corresponding bit in the neigborhood of (i, j) in the image. */
void applyMorphAnd(binaryImage3d &input, 
				   binaryImage3d &morphkernel, 
				   binaryImage3d &mask, 
				   binaryImage3d &output) {
  ulong irows = input.rows();
  ulong icols = input.cols();
  ulong ilays = input.layers();

  float dx = input.getDX();
  float dy = input.getDY();
  float dz = input.getDZ();

  output.setDX(dx);
  output.setDY(dy);
  output.setDZ(dz);

  output.setRows(irows);
  output.setCols(icols);
  output.setLayers(ilays);
  output.resize();

  ulong mrows = morphkernel.rows();
  ulong mcols = morphkernel.cols();
  ulong mlays = morphkernel.layers();
  ulong mrowsOver2 = mrows/2;
  ulong mcolsOver2 = mcols/2;
  ulong mlaysOver2 = mlays/2;
	
  for (int l = 0; l < static_cast<int>(ilays); l++) {
    int ls = l - mlaysOver2;
	for (int i = 0; i < static_cast<int>(irows); i++) {
      int is = i - mrowsOver2;
      for (int j = 0; j < static_cast<int>(icols); j++) {
        int js = j - mcolsOver2;

        bool val = 1;
        for (int ll = max(0, -ls); ll < min(static_cast<int>(mlays), static_cast<int>(ilays) - ls); ll++) {
          for (int ii = max(0, -is); ii < min(static_cast<int>(mrows), static_cast<int>(irows) - is); ii++) {
            for (int jj = max(0, -js); jj < min(static_cast<int>(mcols), static_cast<int>(icols) - js); jj++) {
              if (mask.getVoxel(ii, jj, ll)) {
                val &= (input.getVoxel(is+ii, js+jj, ls+ll) == morphkernel.getVoxel(ii, jj, ll));
              }
            }
          }
        }
        output.setVoxel(i, j, l, val);	
      }
	}
  }
  return;
}


void applyMorphAnd(grayImage3d<GLubyte> &input, 
				   binaryImage3d &morphkernel, 
				   binaryImage3d &mask, 
				   grayImage3d<GLubyte> &output) {
  int irows = input.rows();
  int icols = input.cols();
  int ilays = input.layers();

  float dx = input.getDX();
  float dy = input.getDY();
  float dz = input.getDZ();

  output.setDX(dx);
  output.setDY(dy);
  output.setDZ(dz);
  output.setRows(irows);
  output.setCols(icols);
  output.setLayers(ilays);
  output.resize();
	
  int mrows = morphkernel.rows();
  int mcols = morphkernel.cols();
  int mlays = morphkernel.layers();

  int mrowsOver2 = mrows/2;
  int mcolsOver2 = mcols/2;
  int mlaysOver2 = mlays/2;
	
  for (int l = 0; l < ilays; l++) {
    int ls = l - mlaysOver2;
	for (int i = 0; i <irows; i++) {
      int is = i - mrowsOver2;
      for (int j = 0; j < icols; j++) {
        int js = j - mcolsOver2;

        bool val = 1;
        for (int ll = max(0, -ls); ll < min(mlays, ilays - ls); ll++) {
          for (int ii = max(0, -is); ii < min(mrows, irows - is); ii++) {
            for (int jj = max(0, -js); jj < min(mcols, icols - js); jj++) {
              if (mask.getVoxel(ii, jj, ll)) {
                val &= (input.getVoxel(is+ii, js+jj, ls+ll) == morphkernel.getVoxel(ii, jj, ll));
              }
            }
          }
        }
        if (val > 0) {
          output.setVoxel(i, j, l, 0xff);
        } else {
          output.setVoxel(i, j, l, 0x00);
        }
      }
	}
  }
  return;
}


/* applyMorphOr applies a morphologica operation using the given mask and kernel. The suffix "Or"
 * indicates that the pixel at (i, j) is set if any bit in the kernel, within the support of the mask,
 * matches a corresponding bit in the neigbhorhood of (i, j) in the image. */
void applyMorphOr(binaryImage3d &input, 
				  binaryImage3d &morphkernel, 
				  binaryImage3d &mask, 
				  binaryImage3d &output) {
  int irows = input.rows();
  int icols = input.cols();
  int ilays = input.layers();

  float dx = input.getDX();
  float dy = input.getDY();
  float dz = input.getDZ();

  output.setDX(dx);
  output.setDY(dy);
  output.setDZ(dz);
	
  output.setRows(irows);
  output.setCols(icols);
  output.setLayers(ilays);
  output.resize();
	
  int mrows = morphkernel.rows();
  int mcols = morphkernel.cols();
  int mlays = morphkernel.layers();

  int mrowsOver2 = mrows/2;
  int mcolsOver2 = mcols/2;
  int mlaysOver2 = mlays/2;
	

  for (int l = 0; l < ilays; l++) {
    int ls = l - mlaysOver2;
	for (int i = 0; i < irows; i++) {
      int is = i - mrowsOver2;
      for (int j = 0; j < icols; j++) {
        int js = j - mcolsOver2;

        bool val = 0;
        for (int ll = max(0, -ls); ll < min(mlays, ilays - ls); ll++) {
          for (int ii = max(0, -is); ii < min(mrows, irows - is); ii++) {
            for (int jj = max(0, -js); jj < min(mcols, icols - js); jj++) {
              if (mask.getVoxel(ii, jj, ll)) {
                val |= input.getVoxel(is + ii, js + jj, ls + ll) == morphkernel.getVoxel(ii, jj, ll);
              }
            }
          }
        }
        output.setVoxel(i, j, l, val);
      }
    }
  }
  return;
}

void applyMorphOr(grayImage3d<GLubyte> &input, 
				  binaryImage3d &morphkernel, 
				  binaryImage3d &mask, 
				  grayImage3d<GLubyte> &output) {
  int irows = input.rows();
  int icols = input.cols();
  int ilays = input.layers();

  float dx = input.getDX();
  float dy = input.getDY();
  float dz = input.getDZ();

  output.setDX(dx);
  output.setDY(dy);
  output.setDZ(dz);

  output.setRows(irows);
  output.setCols(icols);
  output.setLayers(ilays);
  output.resize();
	
  int mrows = morphkernel.rows();
  int mcols = morphkernel.cols();
  int mlays = morphkernel.layers();

  int mrowsOver2 = mrows/2;
  int mcolsOver2 = mcols/2;
  int mlaysOver2 = mlays/2;
	

  for (int l = 0; l < ilays; l++) {
    int ls = l - mlaysOver2;
	for (int i = 0; i < irows; i++) {
      int is = i - mrowsOver2;
      for (int j = 0; j < icols; j++) {
        int js = j - mcolsOver2;

        bool val = 0;
        for (int ll = max(0, -ls); ll < min(mlays, ilays - ls); ll++) {
          for (int ii = max(0, -is); ii < min(mrows, irows - is); ii++) {
            for (int jj = max(0, -js); jj < min(mcols, icols - js); jj++) {
              if (mask.getVoxel(ii, jj, ll)) {
                val |= input.getVoxel(is + ii, js + jj, ls + ll) == morphkernel.getVoxel(ii, jj, ll);
              }
            }
          }
        }
        if (val > 0) {
          output.setVoxel(i, j, l, 0xff);
        } else {
          output.setVoxel(i, j, l, 0x00);
        }
      }
    }
  }
  return;
}


void erosionOp(binaryImage3d &input, binaryImage3d &output) {
  applyMorphAnd(input, smallCube, smallCubeMask, output);
}

void erosionOp(grayImage3d<GLubyte> &input, grayImage3d<GLubyte> &output) {
  applyMorphAnd(input, smallCube, smallCubeMask, output);
}

void dialationOp(binaryImage3d &input, binaryImage3d &output) {
  applyMorphOr(input, mediumCube, mediumCubeMask, output);
}

void dialationOp(grayImage3d<GLubyte> &input, grayImage3d<GLubyte> &output) {
  applyMorphOr(input, mediumCube, mediumCubeMask, output);
}

void closureOp() {
  applyMorphOr(afterThreshold, mediumCube, mediumCubeMask, processedImage[0]);
  applyMorphAnd(processedImage[0], mediumCube, mediumCubeMask, processedImage[1]);
}
