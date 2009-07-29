// 
//  main.cc
//  nuclearSurface
//  
//  Created by Robert Snapp on 2007-04-13.
//  Copyright 2007 Robert R. Snapp. All rights reserved.
// 

// Program nuclearSurface performs a three-dimensional analysis of the shape of the largest cell
// nucleus that is contained in a specified Zeiss laser-scanning-microscope (lsm) image file.
// During processing, all voxels with intensities greater than or equal to a prescribed threshold
// value are activated. A morphological closure operation is then performed, in order to fill in
// pores. Using a 26-nearest-neighbor topology, connected groups of voxels are identified and labeled.
// Finally, the convex hull of the largest component is constructed, the ratio of the volume of the
// convex hull to that of the largest component is computed.

#include "lsm.h"
#include "polyhedron.h"
#include "convexhull.h"
#include "imageComponentArray.h"
#include "lsmVoxelImage.h"
#include "ntuple.h"
#include "nuclearSurface.h"
#include "binaryImage3d.h"
#include "grayImage3d.h"
#include "env.h"

#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <climits>

using namespace std;

// Global Variables.
int   const band_def         =  0;
int   const threshold_def     = 15;
int   const verbosity_def     =  0;
char  const FILE_NAME_DELIMS[]    = " \t\r\n,;";

// Define the global variables:
int  band;
bool graphics = false;
string fileName;
string directory;
// string path;
int threshold;
int verbosity;

// Images
lsmVoxelImage    lsmData;					  // The original (base) voxel image
binaryImage3d afterThreshold; 				  // Threshold map of the base image.
vector<binaryImage3d>	 processedImage(2);   // Image buffers for the closure operation.
// grayImage3d<GLubyte> afterThreshold;
// std::vector<grayImage3d<GLubyte> > processedImage(2);

ImageComponentArray comp;					  // Stores labels that identify each connected component.

ConvexHull<int> ch;				// The convex hull of the largest component.

Subvolume subvol = {0, INT_MAX, 0, INT_MAX, 0, INT_MAX};
float dx, dy, dz;

vector<string> fileRoster;                     // The contents of an optional input file that contains a list of
											   // file names.
vector<string>::iterator fileRosterPos;
int fileRosterIndex = 0;

// Output file indicates that path where the results of the convex analysis are printed. If the string
// outputFilename is empty, then results will appear on stdout.
string outputFilename;          
std::ofstream output;

// Function eraseDirectories deletes the names of any leading directories that appear in the path name
// contained in the string path. Thus if path contains the string "/usr/foo/bar.h" then after invoking
// eraseDirectories(path), it will contain "bar.h". This function assumes that the directory separation
// character is a forward slash ('/'). The return value is a reference to the modified path.
string& eraseDirectories(string &path) {
 	size_t last_slash = path.find_last_of("/");  // Remove leading directories
   	if (last_slash != string::npos)
   		(void) path.erase(0, last_slash + 1);
   	return path;
}


// Function intializeGlobals assigns default values for all global variables. It should be invoked
// before the command line arguments are parsed.
void initializeGlobals() {
  band      = band_def;
  threshold = threshold_def;
  verbosity = verbosity_def;
  graphics = false;
}

// applyThreshold performs a threshold operation to each voxel in the indicated channel and subvolume of the
// indicated lsmVoxelImage. The corresponding bit in the binaryImage3d afterThreshold is set to 1 if the voxel
// value exceeds the indicated threshold. The dimensions of afterThreshold are adjusted to correspond to the
// indicated subvolume. 
void applyThreshold(lsmVoxelImage &lsmData, 
                    int channel, 
                    Subvolume &sv, // This structure defines the rectangular volume of interest.
                    int threshold, 
                    binaryImage3d &outImage) {

  if (verbosity > 1) {
    cerr << "Calling applyThreshold with threshold = " 
		 << threshold 
		 << ", channel = " 
		 << channel << "... " 
		 << endl;
  }

  int nRows = sv.height;
  int nCols = sv.width;
  int nLays = sv.depth;

  float dx = lsmData.getDX();  // retrieve the voxel dimensions, in meters.
  float dy = lsmData.getDY();
  float dz = lsmData.getDZ();

  outImage.setDX(dx);  // store the lsm voxel dimensions into the voxel image.
  outImage.setDY(dy);
  outImage.setDZ(dz);

  outImage.setRows(nRows);
  outImage.setCols(nCols);
  outImage.setLayers(nLays);
  outImage.resize();

  for(int k = 0; k < nLays; k++) {
	int kk = sv.layer + k;

	for(int i = 0; i < nRows; i++) {
	  int ii = sv.row + i;

	  for(int j = 0; j < nCols; j++) {
		int jj = sv.col + j;
		bool val = lsmData.getVoxel(ii, jj, kk, channel) >= threshold;
        outImage.setVoxel(i, j, k, val);
      }
    }
  }
  if (verbosity > 1) {
    cerr << " done [" << nRows << "*" << nCols << "*" << nLays << "]" << endl;
  }
}

void applyThreshold(lsmVoxelImage &lsmData, 
                    int channel, 
                    Subvolume &sv, // This structure defines the rectangular volume of interest.
                    int threshold, 
                    grayImage3d<GLubyte> &outImage) {

  if (verbosity > 1) {
    cerr << "Calling applyThreshold with threshold = " 
		 << threshold 
		 << ", channel = " 
		 << channel << "... " 
		 << endl;
  }

  int nRows = sv.height;
  int nCols = sv.width;
  int nLays = sv.depth;

  float dx = lsmData.getDX();
  float dy = lsmData.getDY();
  float dz = lsmData.getDZ();

  outImage.setDX(dx);
  outImage.setDY(dy);
  outImage.setDZ(dz);

  outImage.setRows(nRows);
  outImage.setCols(nCols);
  outImage.setLayers(nLays);
  outImage.resize();

  for(int k = 0; k < nLays; k++) {
	int kk = sv.layer + k;

	for(int i = 0; i < nRows; i++) {
	  int ii = sv.row + i;

	  for(int j = 0; j < nCols; j++) {
		int jj = sv.col + j;
		GLubyte val = (lsmData.getVoxel(ii, jj, kk, channel) >= threshold ? 0xff : 0x00);
        outImage.setVoxel(i, j, k, val);
      }
    }
  }
  if (verbosity > 1) {
    cerr << " done [" << nRows << "*" << nCols << "*" << nLays << "]" << endl;
  }
}

// Function parseToken parses each string token on a fileRoster line.
void parseToken(string const &token, int field) {
  int const FILENAME    = 0;
  int const COL         = 1;
  int const ROW         = 2;
  int const LAYER       = 3;
  int const WIDTH       = 4;
  int const HEIGHT      = 5;
  int const DEPTH       = 6;
  int const THRESHOLD   = 7;

  switch(field) {
  case FILENAME:
    fileName = token;
    break;
  case COL:
    subvol.col = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "subvol.col set to " << subvol.col << endl;
    }
    break;
  case ROW:
    subvol.row = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "subvol.row set to " << subvol.row << endl;
    }
    break;
  case LAYER:
    subvol.layer = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "subvol.layer set to " << subvol.layer << endl;
    }
    break;
  case WIDTH:
    subvol.width = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "subvol.width set to " << subvol.width << endl;
    }
    break;
  case HEIGHT:
    subvol.height = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "subvol.height set to " << subvol.height << endl;
    }
    break;
  case DEPTH:
    subvol.depth = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "subvol.depth set to " << subvol.depth << endl;
    }
    break;
    
  case THRESHOLD:
    threshold = atoi(token.c_str());
    if (verbosity > 1) {
      cerr << "threshold set to " << threshold << endl;
    }
    break;
  default:
    cerr << "Can't parse token " << token
              << " which appears in field " << field
              << endl;
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
      p1 = line.find_first_of(FILE_NAME_DELIMS, p0); // now p1 points to the end
      string token = line.substr(p0, p1);
      parseToken(token, f++);
    }
  }
  return f > 0;
}
    
// Function loadAndPrepareLsmFile reads the indicated lsm file and completes the initial
// (noninteractive) processing.  The current threshold is
// applied to the current layer and the results are stored in the binary image afterThreshold.
bool loadAndPrepareLsmFile(const char* filename, Subvolume &sv) {
  if (verbosity > 1) {
    cerr << "Calling loadAndPrepareLsmFile(" << filename << ", ...)" << endl;
  }

  string path = directory;
  path.append(filename); 

  lsmData.clear();
  bool status = lsmData.read(path);
  if (! status) {
    if (verbosity > 1) {
      cerr << "loadAndPrepareLsmFile(" << filename << "): "
           << "could not access file " << path << endl;
    }
    return false;
  }

  int nRows    = lsmData.rows();
  int nCols    = lsmData.cols();
  int nLayers  = lsmData.layers();
  int nBands   = lsmData.bands();
  // subvolume dimensions that are initially zero should be adjusted
  // to match those of the image.

  if (sv.width == 0 || sv.width + sv.col > nCols) {
    sv.col   = 0;
    sv.width = nCols;
    if (verbosity > 1) {
      cerr << "subvol.width set to " << subvol.width << endl;
    }
  }

  if (sv.height == 0 || sv.height + sv.row > nRows) {
    sv.row    = 0;
    sv.height = nRows;
    if (verbosity > 1) {
      cerr << "subvol.height set to " << subvol.height << endl;
    }
  }

  if (sv.depth == 0 || sv.depth + sv.layer > nLayers) {
    sv.layer = 0;
    sv.depth = nLayers;
    if (verbosity > 1) {
      cerr << "subvol.depth set to " << subvol.depth << endl;
    }
  }

  dx = lsmData.getDx()*microns_per_meter;
  dy = lsmData.getDy()*microns_per_meter;
  dz = lsmData.getDz()*microns_per_meter;

  if (band >= nBands) {
    cout << "Requested band " << band << " does not exist."
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

  applyThreshold(lsmData, band, sv, threshold, afterThreshold);
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
	cerr << "Reading filenames from file " << path << ": ";
  }

  in.open(path.c_str(), ifstream::in);
  if (! in) {
    cerr << endl << "Could not open file " << path << std::endl;
    perror(path.c_str());
    return -1;
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
	cerr << endl 
		 << "(" << count << " lines read from " << path << ")" 
		 << endl;
  }
  return count;
}

// Function analyzeLsmFile completes the batch mode image processing, by first performing
// a closure operation on the image, then labeling the components, computing the convex
// hull of the largest component, and measuring its area and perimeter. The results are
// printed to the indicated output stream.
void analyzeLsmFile(ostream &os) {
  ulong innerHullVolume;  // convex hull measurements are in geometric units.
  ulong innerHullSurface;
  ulong outerHullVolume;
  ulong outerHullSurface;
  double innerRatio;
  double outerRatio;
  
  ulong clusterVolume;
  double dx, dy, dz;
  double cubicMicronsPerVoxel;

  comp.clear();  // Reset the global component array, comp.

  if (verbosity > 5) {
    cout << "Called analyzeLsmFile(...)" << endl;
  }

  // Apply a dialation-erosion sequence to close any isolated dark (off) pixels in the
  // thresholded image (afterThreshold). This generates processedImage[1].
  closureOp();

  // Identify the connected components in the processed image, and label the largest with
  // a unique label. The total number of distinct (isolated) components is returned.
  // Note that an 26-fold topology is used to define the neighborhood of each interior
  // pixel.

  if (verbosity > 5) {
       cout << "\t Calling  comp.computeComponents(processedImage[1])" << endl;
  }
  comp.computeComponents(processedImage[1]);

  dx = comp.getDX();  // get voxel dimensions in meters.
  dy = comp.getDY();
  dz = comp.getDZ();
  cubicMicronsPerVoxel = dx*dy*dz*pow(microns_per_meter, 3);
  

  if (verbosity > 1) {
    cout << "Component scaling factors (in microns), (dx, dy, dz) = (" 
         << dx*microns_per_meter << ", "
         << dy*microns_per_meter << ", "
         << dz*microns_per_meter << ")" << endl;
  }

  
   clusterVolume = comp.clusterSize(0);
  if (verbosity > 1) {
	cout << "comp.clusterSize(0) = " << clusterVolume << " voxels = "
		 << clusterVolume*cubicMicronsPerVoxel << " cubic microns." << endl;
  }

  comp.computeConvexHull(0, ch, VOXEL_CENTERS);
  innerHullVolume = ch.volume();
  innerHullSurface = ch.surfaceArea();
  innerRatio =  static_cast<double>(innerHullVolume) / static_cast<double>(clusterVolume);

  if (verbosity > 1) {
	cout << "inner Hull volume = " 
		 << innerHullVolume << " voxels = " 
		 << innerHullVolume*cubicMicronsPerVoxel << " cubic microns." 
		 << endl;
	cout << "inner hull contains" << ch.vertices() << " vertices, "
		 << ch.edges() << " edges, "
		 << ch.faces() << " faces." << endl;
	cout << "inner hull ration = " 
		 << innerHullVolume << " / " << clusterVolume << " = "
		 << innerRatio << endl;
	if (verbosity > 50) {
	  ch.describe(cout);
	}
  }
  
  //ch.clear();
  comp.computeConvexHull(0, ch, VOXEL_CORNERS);
  outerHullVolume = ch.volume();
  outerHullSurface = ch.surfaceArea();
  outerRatio =  static_cast<double>(outerHullVolume) / static_cast<double>(clusterVolume);

  if (verbosity > 1) {
	cout << "outer Hull volume = " 
		 << outerHullVolume << " voxels = " 
		 << outerHullVolume*cubicMicronsPerVoxel << " cubic microns." 
		 << endl;
  	cout << "outer hull contains" << ch.vertices() << " vertices, "
		 << ch.edges() << " edges, "
		 << ch.faces() << " faces." << endl;
	cout << "outer hull ration = " 
		 << outerHullVolume << " / " << clusterVolume << " = "
		 << outerRatio << endl;
	if (verbosity > 50) {
	  ch.describe(cout);
	}
  }
   
 

  // TODO: Add graphics option to view convex hull in 3D.
	
 
  // Print the results of the convexity analysis of the largest component to the indicated
  // output stream.
  if (os) {
	os << fileName << ", " 
	   << subvol.col << ", "
	   << subvol.row << ", "
	   << subvol.layer << ", "
	   << subvol.width << ", "
	   << subvol.height << ", "
	   << subvol.depth << ", "
	   << band << ", "
	   << threshold << ", "
	   << clusterVolume << ", "
	   << innerHullVolume << ", "
	   << innerHullSurface << ", "
	   << outerHullVolume << ", "
	   << outerHullSurface << ", "
	   << innerRatio << ", "
	   << outerRatio
	   << endl;
  }
}


int main (int argc, char** argv) {
  initializeGlobals();
  initializeGlut(&argc, argv);

  initKernels(); // initialize the kernels in morphology.cc

  int nOpts = parseCommandLine(argc, argv);

  if (! outputFilename.empty()) {
    if(! openOFStream(outputFilename, output)) {
      cerr << argv[0]
           << ": ERROR: can't write to file " << outputFilename
           << endl;
      abort();
	} else {
	  output << "fileName, col, row, layer, width, height, depth, "
			 << "band, threshold, clusterVolume, inner-hull volume, inner-hull surface area, "
			 << "outer-hull volume, outer-hull surface area, "
			 << "inner ratio (inner-hull vol/cluster vol), outer ratio (outer-hull vol/cluster vol)"
			 << endl;
	}
  }
  
  if (verbosity > 1) {
    cout << "band = " << band << ", "
         << "directory = " << directory << ", "
         << "row   = " << subvol.row << ", "
         << "height   = " << subvol.height << ", "
         << "col   = " << subvol.col << ", "
         << "width   = " << subvol.width << ", "
         << "layer   = " << subvol.layer << ", "
         << "depth   = " << subvol.depth << ", "
         << "threshold = " << threshold << ", "
         << "verbosity = " << verbosity << ", "
         << "graphics (mode)  = " << graphics << ". "
         << endl;
  }

  if (argc < nOpts + 1) {
    //Error: No path found
    string command = argv[0];          // Identify the command name
    (void) eraseDirectories(command);  // Remove leading directories (if present).
   		 
    describeUsage(command.c_str());
    return EXIT_FAILURE;
  }
  
  // string path = argv[optind];
  // fileName = eraseDirectories(path);
  cout << argv[nOpts] << endl;
  fileName.assign(directory);
  fileName.append(argv[nOpts]);

  if (verbosity > 1) {
    cout << "nOpts = " << nOpts << ", path = " << fileName << endl; 
  }

  if (isLSM_File(fileName.c_str())) {
    loadAndPrepareLsmFile(fileName.c_str(), subvol);
	if (! graphics) {
     if (output.is_open()) {
        analyzeLsmFile(output);
      } else {
        analyzeLsmFile(cout);
      }
	} else {
	  runGlutLoop();
	}
  } else {
	// The last argument should indicate a file that contains a list of lsm file paths.
    // Open the file, and transcribe the lines to the vector<string> fileRoster.

	int fileCount = readLinesFromFile(fileName.c_str(), fileRoster);
    if (fileCount < 1) {
      cerr << "Cannot process the roster file " << fileName.c_str() << endl;
      return EXIT_FAILURE;
    }
    
    if (verbosity > 1) {
      cerr << "Begining of main loop." << endl;
    }

	if (! graphics) {
	  for(fileRosterPos = fileRoster.begin(); fileRosterPos != fileRoster.end(); ++fileRosterPos) {
		parseFileLine(*fileRosterPos);
		loadAndPrepareLsmFile(fileName.c_str(), subvol);
		if (output.is_open()) {
		  analyzeLsmFile(output);
		} else {
		  analyzeLsmFile(cout);
		}
	  }
	} else {
	  // if graphics mode is enabled then just queue up the first file
      // and start glut interactive loop in multiple file mode.
      fileRosterPos = fileRoster.begin();

	  while (! parseFileLine(*fileRosterPos)) {
		if (fileRosterPos != fileRoster.end()) {
		  ++fileRosterPos;
		} else {
		  cerr << "Could not find a valid line in the file roster to parse." << endl;
		  exit(1);
		}
	  }
	
	  loadAndPrepareLsmFile(fileName.c_str(), subvol);
	  runGlutLoop();
    }
  }
 

  // Identify the largest connected component in the lsm image.
  // TODO: Add z-min and z-max arguments on the command line.
  
  return 0;
}
