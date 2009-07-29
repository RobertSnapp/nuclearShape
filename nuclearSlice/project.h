#ifndef __PROJECT_H__
#define __PROJECT_H__

/****************************************************************************
 * 
 * project.h
 * Project: nuclearSlice
 *
 * Created by Robert Snapp on 2007-07-17
 *
 * Contains data structure definitions, declarations of public functions, and 
 * declarations of global variables required for the nuclearSlice program.
 * 
 ******************************************************************************/

#include "convexHull2d.h"
#include "binaryImage.h"
#include "lsmVoxelImage.h"
#include "env.h"   // for path to glut.h 
#include "grayImage.h"

#include <utility>  // for  pair stl template.
#include <string>
#include <vector>
#include <ostream>
#include <fstream>
#include <iostream>

typedef unsigned long ulong;

struct Rect {
  int row;
  int col;
  int height;
  int width;
};

struct Cluster {
	int size;
	double center[2];
	double lam1;
	double lam2;
	double ev1[2];
	double ev2[2]; 
};

// Defined in main.cc:
void applyThreshold(Rect &region);
void voxelProjection(lsmVoxelImage &v, Rect &region, grayImage<GLint> &g);
void rescaleImage(grayImage<GLint> &x, grayImage<GLubyte> &y);
void applyThresholdToLsmData(int t, Rect &region);
void applyThresholdToScaledImage(int t);
void parseToken(std::string const &token, int field);
bool parseFileLine(std::string const &line);
bool loadAndPrepareLsmFile(const char* filename, Rect &subregion);

// Defined in cmdline.cc:
void describeUsage(char const *command);
int parseCommandLine(int argc, char** argv);

// Defined in components.cc:
int collectComponents();
void computeComponentHull(int label, ConvexHull2d<double> &h);
int countComponents(grayImage<GLint> &labels, std::vector< std::pair<int,int> > &components);
bool greaterFrequency(std::pair<int,int> lhs, std::pair<int,int> rhs);
void mergeComponents(grayImage<GLint> &a, ulong r, ulong c, int from, int to);
int labelComponents(binaryImage &input, grayImage<GLint> &labels);
int labelComponents(grayImage<GLubyte> &input, grayImage<GLint> &labels);
void makeColorLabels(int nLabels);
void computeEigenvalues(int label, Cluster &blob);

// Defined in graphicalInterface.cc:
void initializeGlut(int* argc, char** argv);
void runGlutLoop();

// Defined in morphology.cc
void erosionOp(binaryImage &input, binaryImage &output);
void erosionOp(grayImage<GLubyte> &input, grayImage<GLubyte> &output);
void dialationOp(binaryImage &input, binaryImage &output);
void dialationOp(grayImage<GLubyte> &input, grayImage<GLubyte> &output);
void closureOp();
void initKernels();

/* Global data */
extern bool graphics;
extern int verbosity;
extern int threshold;
extern int layer;
extern int band;
extern std::string directory;

extern lsmVoxelImage lsmData;                       // The original voxel image
extern grayImage<GLint> projectedImage;             // The projection of the voxels (computed in projectionMode)
extern grayImage<GLubyte> scaledImage;              // The projectedImage scaled to an 8-bit depth.

// An artifact occurs when trying to render large (e.g. 1024 x 1024) binary images
// on my Macintosh Powerbook. The left half of the image is repeated in the right
// half of the view. Consequently, the application can be built so that the threshold and
// other binary images are implemented as 8-bit graylevel images using 
// 0x00 and 0xff as the two states. Very wasteful, but the geometry comes
// out correct.

#ifdef EIGHTBIT_MORPHOLOGY_BUFFERS
extern grayImage<GLubyte> afterThreshold;                  // Either the threshold of an lsmData layer, or of the scaledImage.
extern std::vector<grayImage<GLubyte> > processedImage;      // Image buffers for the closure operation.
#else // BINARY_MORPHOLOGY_BUFFERS are assumed by default. Works on the MacBook Pro, but not on the Powerbook G4.
extern binaryImage afterThreshold;                  // Either the threshold of an lsmData layer, or of the scaledImage.
extern std::vector<binaryImage> processedImage;      // Image buffers for the closure operation.
#endif

extern grayImage<GLint> componentLabels;            // Buffer to hold the labels of each component.
extern grayImage<GLubyte> colorLabels;              // Buffer to hold colors of each component.
extern std::vector<std::pair<int, int> > labelCount;          // The size of each component, e.g., the frequency of each label.
extern ConvexHull2d<double> ch;                     // The convex hull of the largest component.
extern Cluster blob;

extern Rect subregion;                              // The rectangular region of interest in all of the above images.

extern std::string fileName;
extern Rect subregion;
extern bool projectionMode;
extern std::vector<std::string> fileRoster;
extern std::vector<std::string>::iterator fileRosterPos;

// Output file indicates that path where the results of the convex analysis are printed. If the string
// outputFilename is empty, then results will appear on stdout.
extern std::string outputFilename;          
extern std::ofstream output;

#endif  //  __PROJECT_H__
