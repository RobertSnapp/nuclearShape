// 
//  main.cc
//  nuclearSlice
//  
//  Created by Robert Snapp on 2007-04-13.
//  Copyright 2007 Robert R. Snapp. All rights reserved.
// 

// Program nuclearSlice automatically describes the shape of cell nucleii, and other particles, that are imaged
// in one or more lsm files. In graphical mode (enabled with the command line option -G) a user can interactively
// adjust a variety of viewing parameters (image layer, image band, threshold) before applying the principle
// image processing operation in which each group of adjacent pixels above threshold is identified as a distinct
// component, and the convex hull of the largest component is then constructed, and its area computed. In batch
// mode the viewing parameters can be specified on the command line, and the results are printed on the standard
// output stream.
//
// Project organization:
//
// project.h                contains the declarations of functions required by other modules, and declarations of
//                          global variables.
// graphicalInterface.cc    defines functions that support the interactive graphical mode. Currently, glut and openGL
//                          are used.
// components.cc            defines functions for grouping adjacent pixels above threshold into distinct components.
// cmdline.cc               defines a function to parse the command line
// morphology.cc            defines functions for simple image morphology operations: dialation, erosion, closure.
// main.cc                  defines the main program procedure, as well as auxillary image processing operations.
//

#include "project.h"
#include "env.h"
#include "lsm.h"
#include "polygon.h"
#include "convexHull2d.h"
#include "binaryImage.h"
#include "grayImage.h"
#include "lsmVoxelImage.h"
#include "ntuple.h"

#ifdef PNG_IMAGE
#include "pngImage.h"
#endif

#include <algorithm>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include <cstdlib>
#include <utility>  // for pair
#include <string>

using namespace std;

void eraseDirectories(string &path);

// Global Variables.
int   const band_def         =  0;
int   const layer_def         =  0;
int   const threshold_def     = 40;
int   const verbosity_def     =  0;
float const microns_per_meter = 1.0E06;
char  const FILE_NAME_DELIMS[]    = " \t\n,;";

bool graphics;
bool projectionMode;
int  verbosity;
int  threshold;
int  layer;
int  band;
string fileName("filename");

// directory can be set by the command line option -d. Its value will then be prepended to the path
// of each lsm file.
string directory;   
 

// Images
lsmVoxelImage    lsmData;                    // The original voxel image
grayImage<GLint> projectedImage;             // The projection of the voxels (computed in projectionMode)
grayImage<GLubyte> scaledImage;              // The projectedImage scaled to an 8-bit depth.

// Because of an artifact that occurs when rendering a binary image on my macintosh powerbook (the left
// half of the image is also drawn over the right half of the view), binary images will be represented
// using a grayscale format, with 0x00 for black, and 0xff for white.

#ifdef EIGHTBIT_MORPHOLOGY_BUFFERS
grayImage<GLubyte> afterThreshold;                  // Either the threshold of an lsmData layer, or of the scaledImage.
std::vector<grayImage<GLubyte> > processedImage;      // Image buffers for the closure operation.
#else // BINARY_MORPHOLOGY_BUFFERS are assumed by default. Works on the MacBook Pro, but not on the Powerbook G4.
binaryImage afterThreshold;                  // Either the threshold of an lsmData layer, or of the scaledImage.
std::vector<binaryImage> processedImage;      // Image buffers for the closure operation.
#endif


grayImage<GLint> componentLabels;            // Buffer to hold the labels of each component.
grayImage<GLubyte> colorLabels;              // Buffer to hold colors of each component.
vector<pair<int, int> > labelCount;          // The size of each component, e.g., the frequency of each label.
Cluster blob;                                // Data structure defined in project.h for characterizing the
                                             // elliptical eccentricity of the largest blob.
ConvexHull2d<double> ch;                     // The convex hull of the largest component.

Rect subregion;                              // The rectangular region of interest in all of the above images.
float dx, dy, dz;

vector<string> fileRoster;                   // The contents of an optional input file that contains a list of
                                             // file names.
vector<string>::iterator fileRosterPos;
int fileRosterIndex = 0;

// Output file indicates that path where the results of the convex analysis are printed. If the string
// outputFilename is empty, then results will appear on stdout.
string outputFilename;          
ofstream output;


// Function eraseDirectories deletes the names of any leading directories that appear in the path name
// contained in the string path. Thus if path contains the string "/usr/foo/bar.h" then after invoking
// eraseDirectories(path), it will contain "bar.h". This function assumes that the directory separation
// character is a forward slash ('/'). The return value is a reference to the modified path.
void eraseDirectories(string &path) {
  size_t last_slash = path.find_last_of("/");  // Remove leading directories
  if (last_slash != string::npos)
    (void) path.erase(0, last_slash + 1);
  return;
}


// Function intializeGlobals assigns default values for all global variables. It should be invoked
// before the command line arguments are parsed.
void initializeGlobals() {
  band      = band_def;
  layer     = layer_def;
  threshold = threshold_def;
  verbosity = verbosity_def;
  graphics = false;
  projectionMode = false;
  subregion.row = 0;
  subregion.col = 0;
  subregion.width = 0;
  subregion.height = 0;
  
#ifdef EIGHTBIT_MORPHOLOGY_BUFFERS
  processedImage.push_back(grayImage<GLubyte>());
  processedImage.push_back(grayImage<GLubyte>());
#else // BINARY_MORPHOLOGY_BUFFERS are assumed by default. Works on the MacBook Pro, but not on the Powerbook G4.
  processedImage.push_back(binaryImage());
  processedImage.push_back(binaryImage());
#endif
}

// Function parseToken parses each string token on a fileRoster line.
void parseToken(string const &token, int field) {
  int const FILENAME = 0;
  int const COL = 1;
  int const ROW = 2;
  int const WIDTH = 3;
  int const HEIGHT = 4;
  int const THRESHOLD = 5;

  switch(field) {
  case FILENAME:
	fileName.assign(directory);
    fileName.append(token);
    break;
  case COL:
    subregion.col = atoi(token.c_str());
    break;
  case ROW:
    subregion.row = atoi(token.c_str());
    break;
  case WIDTH:
    subregion.width = atoi(token.c_str());
    break;
  case HEIGHT:
    subregion.height = atoi(token.c_str());
    break;
  case THRESHOLD:
    threshold = atoi(token.c_str());
    break;
  default:
    std::cerr << "Can't parse token " << token
              << " which appears in field " << field
              << std::endl;
    break;
  }
  return;
}

bool parseFileLine(string const &line) {
  // The first word on each line of fileRoster should contain the path of an
  // accessible lsm file.
  int f  = 0; // field index
  size_t p0 = 0; // pointer to token beginning
  size_t p1 = 0; // pointer to token end

  while(p0 != string::npos && p1 != string::npos) {
    p0 = line.find_first_not_of(FILE_NAME_DELIMS, p1);
    if (p0 != string::npos) {
      // p0 points to the beginning of a token.
      p1 = line.find_first_of(FILE_NAME_DELIMS, p0); // now p1 points to the end of the token
      string token = line.substr(p0, p1);
      parseToken(token, f++);
    }
  }
  return f > 0;
}
    
// Function loadAndPrepareLsmFile reads the indicated lsm file and completes the initial
// (noninteractive) processing. In projection mode (projectionMode == true) this involves
// loading the projectedImage buffer with the sum of pixel intensities in each layer of
// the lsm file, scaling the intensities (from 0 to 255) and storing the scaled values in
// scaledImage, and then applying the current threshold to create the binary image
// afterThreshold. In layer mode (profileMode == false), then the current threshold is
// applied to the current layer and the results are stored in the binary image afterThreshold.
// The first argument dir can contain the path to a different directory. Its value is prepended
// to the second argument filename in order to obtain a complete path name. The third argument, sr.
// defines a rectangular subregion of interest with respect to which the lsm data will be clipped.
bool loadAndPrepareLsmFile(const char* filename, Rect &sr) {
  string path = filename;

  lsmData.clear();
  bool status = lsmData.read(path);
  if (! status) {
    if (verbosity > 1) {
      cerr << "loadAndPrepareLsmFile(" << filename << ", ...): "
           << "could not access file " << path << endl;
    }
    return false;
  }

  int nRows    = lsmData.rows();
  int nCols    = lsmData.cols();
  int nLayers  = lsmData.layers();
  int nBands   = lsmData.bands();
  // subregion dimensions that are initially zero should be adjusted
  // to match those of the image.

  if (sr.width == 0 || sr.width + sr.col > nCols) {
    sr.col   = 0;
    sr.width = nCols;
  }

  if (sr.height == 0 || sr.height + sr.row > nRows) {
    sr.row    = 0;
    sr.height = nRows;
  }

  dx = lsmData.getDx()*microns_per_meter;
  dy = lsmData.getDy()*microns_per_meter;
  dz = lsmData.getDz()*microns_per_meter;

  if (layer >= nLayers) {
    cerr << "Requested layer " << layer << " does not exist."
         << endl;
    layer = 0;
  }

  if (band >= nBands) {
    cerr << "Requested band " << band << " does not exist."
         << endl;
    band = 0;
  }

  if (verbosity > 1) {
    cout << "dx = " << dx << " microns, "
         << "dy = " << dy << " microns, "
         << "dz = " << dz << " microns." << endl;
    cout << "rows   = " << nRows    << endl
         << "cols   = " << nCols    << endl
         << "layers = " << nLayers  << endl
         << "bands  = " << nBands   << endl << endl;
  }

  // Compute the threshold map of the current lsm image.
  applyThreshold(sr);

  // Compute the closure map
  closureOp();

  return true;
}


bool openOFStream(string &path, ofstream &os) {
  // An output file is created to record the results
  os.open(path.c_str());
  if (! os) {
    std::perror(path.c_str());
    return false;
  }
	
  return true;
}

// Function readLinesFromFile opens the indicated file, and pushes lines onto the
// back of the indicated string vector. Lines that begin with a comment character
// or newline are ignored. If this task can be completed, the number of retained
// lines is returned. A value of -1 is returned if the file cannot be opened.
int readLinesFromFile(const char* filename, std::vector<string> &lines) {
  string path = filename;
  char buffer[BUFSIZ];
  std::ifstream in;
  int count = 0;

   if (verbosity > 1) {
    cerr << "Reading filenames from file " << path  << ":" << endl;
  }

  in.open(path.c_str(), ifstream::in);
  if (! in) {
    std::cerr << "Could not open file " << path << std::endl;
    std::perror(path.c_str());
    return -1;
  } else if (verbosity > 2) {
	cerr << "Opened file " << path << " for reading." << endl;
  }
	
  while(in.getline(buffer, BUFSIZ)) {
	if (verbosity > 2) {
	  cerr << ".";
	}

    // Skip if this line is blank or is a comment. Comment lines have a leading # character.
    if (buffer[0] != '#' && buffer[0] != '\n') {
      count++;
      lines.push_back(string(buffer));
    }
  }

  in.close();
  if (verbosity > 1) {
    cerr << "(" << count << " filenames read)" << endl;
  }
  return count;
}

// Function analyzeLsmFile completes the batch mode image processing, by first performing
// a closure operation on the image, then labeling the components, computing the convex
// hull of the largest component, and measuring its area and perimeter. The results are
// printed to the indicated output stream.
void analyzeLsmFile(ostream &os) {

  if (verbosity > 5) {
    cerr << "analyzeLsmFile(...)" << endl;
  }

  // Apply a dialation-erosion sequence to close any isolated dark (off) pixels in the
  // thresholded image (afterThreshold).
  closureOp();

  // Identify the connected components in the processed image, and label the largest with
  // a unique label. The total number of distinct (isolated) components is returned.
  // Note that an eight-fold topology is used to define the neighborhood of each interior
  // pixel.
  int nComponents = collectComponents();

  // Print the results of the convexity analysis of the largest component to the indicated
  // output stream.
  os << fileName << ", " 
     << subregion.col << ", "
     << subregion.row << ", "
     << subregion.width << ", "
     << subregion.height << ", ";

  if (projectionMode) {
    os << "P" << lsmData.layers() << ", ";
  } else {
    os << layer << ", ";
  }
  os << band << ", "
     << threshold << ", ";

  if (nComponents > 0) {
	// Compute the convex hull of the largest connected component.
	computeComponentHull(1, ch); 

	if (verbosity > 1) {
	  ch.describe(cerr);
	}
 
	for(size_t i = 0; i < 3 ; ++i) {
	  int nextComponentSize = (i < labelCount.size() ? labelCount[i].second : 0);
	  os << nextComponentSize << ", ";
	}

	os << labelCount[0].second*dx*dy << ", ";  // Physical area of largest blob

	// Compute eignevalues of largest blob
	computeEigenvalues(1, blob);
	double angle = (180./M_PI)*atan2(blob.ev1[1], blob.ev1[0]);

	os << blob.lam1 << ", " << blob.lam2 << ", " << sqrt(blob.lam1/blob.lam2) << ", ";
	os << angle << ", "; // orientation of direction of elongation with respect to positive x-axis.
	os << ch.area() << ", "
	   << ch.perimeter() << ", "
	   << ch.area()/labelCount[0].second
	   << endl;
  } else {
	os << "Warning: no components found." << endl;
  }
}

// Function voxelProjection computes the projection (i.e., the vector sum) of the 8-bit layers in
// the first image (v), and returns that sum in the second (32-bit) image g.
void voxelProjection(lsmVoxelImage &v, Rect &region, grayImage<GLint> &g) {
  g.setRows(region.height);
  g.setCols(region.width);
  g.resize();
  int nLayers = v.layers();

  for(int i = 0, ii = region.row; i < region.height; i++, ii++) {
    for(int j = 0, jj = region.col; j < region.width; j++, jj++) {
      uint32 x = 0;
      for(int k = 0; k < nLayers; k++) {
        x += static_cast<uint32>(v.getVoxel(ii, jj, k, band));
      }
      g.setPixel(i, j, x);
    }
  }
}

// Function rescaleImage converts a 32-bit, grayscale image into an 8-bit gray scale image
// by means of an affine transformation that maps the maximum observed intensity in the first
// to 255 in the second, and the minimul observed intensity in the first, to 0 in the second.
void rescaleImage(grayImage<GLint> &x, grayImage<GLubyte> &y) {
  int nRows = x.rows();
  int nCols = x.cols();
  
  int maxValue = 0;
  int minValue = 0;
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      int z = x.getPixel(i, j);
      maxValue = max(maxValue, z);
      minValue = min(minValue, z);
    }
  }

  y.setRows(nRows);
  y.setCols(nCols);
  y.resize();

  float b = (maxValue > minValue ? 255.0/(maxValue - minValue) : 0);
  float a = - b*minValue;

  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      int z  = x.getPixel(i, j);
      int zp =  static_cast<int>(a + b*z);
      y.setPixel(i, j, zp);
    }
  }

  if (verbosity > 5) {
	if (verbosity > 10) {
	  cerr << "[" << __FILE__ << ":" << __LINE__ << "] ";
	}
    cerr << "rescaleImage(): minValue = " << dec << minValue 
         << ", maxValue = " << maxValue
         << ", a = " << a
         << ", b = " << b
         << endl;
  }

}

// applyThreshold is a utility routine that computes the threshold in a manner that
// is consistent with the current projectionMode. If projectionMode is set, then
// the a projection of the current band of the lsmData is computed, and rescaled.
// Otherwise, the projection is simply applied to the current layers in the lsmData.
void applyThreshold(Rect &region) {
  if (projectionMode) {
	voxelProjection(lsmData, region, projectedImage);
	rescaleImage(projectedImage, scaledImage);
	applyThresholdToScaledImage(threshold);
  } else {
	applyThresholdToLsmData(threshold, region);
  }
}

// Function applyThresholdToLsmData (with two arguments) fills the afterThreshold buffer with a binary image that
// represents the pixels in the current layer and band that exceeds the specified threshold. Note that
// second argument (region) defines the rectangular domain of the thresholded image.
void applyThresholdToLsmData(int t, Rect &region) {
  afterThreshold.setRows(region.height);
  afterThreshold.setCols(region.width);
  afterThreshold.resize();

  for(int i = 0, ii = region.row; i < region.height; i++, ii++) {
    for(int j = 0, jj = region.col; j < region.width; j++, jj++) {

#ifdef EIGHTBIT_MORPHOLOGY_BUFFERS
      GLubyte x = (lsmData.getVoxel(ii, jj, layer, band) > t ? 0xff : 0x00);
#else
      bool x = (lsmData.getVoxel(ii, jj, layer, band) > t);
#endif

      afterThreshold.setPixel(i, j, x);
    }
  }
}

// Function applyThresholdToScaledImage (with two arguments) fills the afterThreshold buffer with a binary image that
// represents the pixels in the scaled projection exceeds the specified threshold. Note that
// second argument (region) defines the rectangular domain of the thresholded image.
void applyThresholdToScaledImage(int t) {
  int nRows = scaledImage.rows();
  int nCols = scaledImage.cols();

  afterThreshold.setRows(nRows);
  afterThreshold.setCols(nCols);
  afterThreshold.resize();

  for(int i = 0; i < nRows; i++) {
    for(int j = 0; j < nCols; j++) {

#ifdef EIGHTBIT_MORPHOLOGY_BUFFERS
      GLubyte x = (scaledImage.getPixel(i, j) > t ? 0xff : 0x00);
#else
      bool x = (scaledImage.getPixel(i, j) > t);
#endif
      afterThreshold.setPixel(i, j, x);
    }
  }
}

int main (int argc, char** argv) {
  initializeGlobals();
  initKernels();    // initialize the kernels in morphology.cc
  initializeGlut(&argc, argv);
  int nOpts = parseCommandLine(argc, argv);
  
  if (! outputFilename.empty()) {
    if(! openOFStream(outputFilename, output)) {
      cerr << argv[0]
           << ": ERROR: can't write to file " << outputFilename
           << endl;
      abort();
    } else {
	  output << "fileName, col, row, width, height, layer (P = projection), "
			 << "band, threshold, c1-size, c2-size, c3-size, c1-area(sq. microns), "
			 << "ev1, ev2, eccentricity, orientation, ch-area, ch-perim, ch-area/c1-size"
			 << endl;
	}

  }

  if (verbosity > 1) {

#ifdef EIGHTBIT_MORPHOLOGY_BUFFERS
	cout << "Using eightbit (unsigned char) morphology buffers." << endl;
#else
	cout << "Using binary morphology buffers." << endl;
#endif

    cout << "band = " << band << ", "
         << "layer   = " << layer << ", "
         << "threshold = " << threshold << ", "
         << "verbosity = " << verbosity << ", "
         << "graphics  = " << graphics << ". "
         << endl;
  }

  // The next argument should be either the path to an lsm file, or the path to a text file that contains a list
  // of lsm filenames.

  if (argc < nOpts + 1) {
    // Usage error: No path was found.
    string command = argv[0];          // Identify the command name
    (void) eraseDirectories(command);  // Remove leading directories (if present).
   	
    cerr << command << ": Usage Error: <filename> argument expected." << endl;

    describeUsage(command.c_str());
    return EXIT_FAILURE;
  }
  cout << argv[nOpts] << endl;
  fileName.assign(directory);
  fileName.append(argv[nOpts]);
  // eraseDirectories(fileName);
  
  if (verbosity > 1) {
    cout << "nOpts = " << nOpts << ", path = " << fileName << endl;
  }

  if (isLSM_File(fileName.c_str())) {
	  loadAndPrepareLsmFile(fileName.c_str(), subregion);
    if (! graphics) {
      if (output.is_open()) {
        analyzeLsmFile(output);
      } else {
        analyzeLsmFile(cout);
      }
    } else {
      // queue up the image for interactive analysis.
      // launch glut interactive loop in single file mode.
      runGlutLoop();
    }
  } else {
    // The last argument should indicate a file that contains a list of lsm file paths.
    // Open the file, and transcribe the lines to the vector<string> fileRoster.
	  int fileCount = readLinesFromFile(fileName.c_str(), fileRoster);
    if (fileCount < 1) {
      cerr << "Cannot process the roster file " << fileName.c_str() << endl;
      exit(1);
    }

    if (! graphics) {
      // If the graphics option "-G" was not specified in the command line argument, then apply
      // the standard sequence of operations to each entry in the file roster, as a batch process,
      // without user intervention.
      for(fileRosterPos = fileRoster.begin(); fileRosterPos != fileRoster.end(); ++fileRosterPos) {
        if (parseFileLine(*fileRosterPos)) {
		  loadAndPrepareLsmFile(fileName.c_str(), subregion);
		  if (output.is_open()) {
			analyzeLsmFile(output);
		  } else {
			analyzeLsmFile(cout);
		  }
		}
      }
    } else {
      // if graphics mode is enabled then just queue up the first file
      // and start glut interactive loop in multiple file mode.
      fileRosterPos = fileRoster.begin();
	  bool status;
	  while (! parseFileLine(*fileRosterPos)) {
		if (fileRosterPos != fileRoster.end()) {
		  ++fileRosterPos;
		} else {
		  cerr << "Could not find a valid line in the file roster to parse." << endl;
		  exit(1);
		}
	  }
	
	  loadAndPrepareLsmFile(fileName.c_str(), subregion);
	  runGlutLoop();
    }
  }
 
  return 0;
}
