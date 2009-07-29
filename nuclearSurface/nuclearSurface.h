/* @(#)nuclearSurface.h
 * Project nuclearSurface
 *
 * Header file for global declarations.
 * Copyright (C) 2007 Robert R. Snapp. All rights reserved.
 */

#ifndef __NUCLEARSURFACE_H__
#define __NUCLEARSURFACE_H__

#include <string>
#include <vector>
#include "binaryImage3d.h"
#include "grayImage3d.h"
#include "imageComponentArray.h"
#include "convexhull.h"
#include "lsm.h"
#include "lsmVoxelImage.h"

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <climits>

class ImageComponentArray;

struct Subvolume {
  int row;
  int height;
  int col;
  int width;
  int layer;
  int depth;
};

// A structure for managing a rectangular region. Used for selection region
// in graphicsInterface.cc, and might be useful elsewhere.
struct Rect {
  int col;
  int row;
  int width;
  int height;
};

// Some useful constants:
const float microns_per_meter = 1.0E06;

// Defined in cmdline.cc
extern int parseCommandLine(int argc, char** argv);
extern void describeUsage(char const *command);

// Defined in graphicalInterface.cc
void initializeGlut(int* argc, char** argv);
void runGlutLoop();

// Defined in main.cc
void applyThreshold(lsmVoxelImage &lsmData, 
                    int channel, 
                    Subvolume &sv, // This structure defines the rectangular volume of interest.
                    int threshold, 
                    binaryImage3d &afterThreshold);
void applyThreshold(lsmVoxelImage &lsmData, 
                    int channel, 
                    Subvolume &sv, // This structure defines the rectangular volume of interest.
                    int threshold, 
                    grayImage3d<GLubyte> &afterThreshold);

void parseToken(std::string const &token, int field);

bool parseFileLine(std::string const &line);
bool loadAndPrepareLsmFile(const char* filename, Subvolume &sv);



// Defined in morphology.cc
void erosionOp(binaryImage3d &input, binaryImage3d &output);
void erosionOp(grayImage3d<GLubyte> &input, grayImage3d<GLubyte> &output);
void dialationOp(binaryImage3d &input, binaryImage3d &output);
void dialationOp(grayImage3d<GLubyte> &input, grayImage3d<GLubyte> &output);
void closureOp();
void initKernels();

extern int band;                      // The image band (or channel) to be analyzed.
extern bool graphics;                 // If true, then convex hull is rendered with OpenGL
extern std::string outputFilename;    // If non-null, then output results are printed in outputFilename
extern std::string directory;         // The directory (path) that contains the lsm file.
extern Subvolume subvol;              // The subvolume of interest: voxels outside the subvolume are ignored.
extern int threshold;                 // Voxels with intensities below this value will be ignored.
extern int verbosity;                 // verbosity controls the volume of diagnostic messages.

// Either the threshold of an lsmData layer, or of the scaledImage.
extern lsmVoxelImage lsmData;
extern binaryImage3d afterThreshold;                  
extern std::vector<binaryImage3d>	 processedImage;      // Image buffers for the closure operation.
// extern grayImage3d<GLubyte> afterThreshold;
// extern std::vector<grayImage3d<GLubyte> > processedImage;

extern ImageComponentArray comp;
extern ConvexHull<int> ch;
extern std::string fileName;
extern std::vector<std::string> fileRoster;
extern std::vector<std::string>::iterator fileRosterPos;

#endif /* __NUCLEARSURFACE_H__ */
