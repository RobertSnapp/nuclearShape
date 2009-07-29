/*
 * components.cc
 * Project nuclearSlice
 *
 * Created by Robert Snapp on 2007-07-16
 *
 * This file contains routines and functions for manipulating sets of
 * adjancent active pixels in a binary image.
 *
 */

#include "project.h"
#include "grayImage.h"
#include "binaryImage.h"
#include <utility> // for pair
#include <vector>

using namespace std;

/* countComponents computes the size of each connected component in the image. The zero-level
 * background is ignored.
 */
int countComponents(grayImage<GLint> &labels, vector< pair<int,int> > &components) {
  if (verbosity > 1) {
     cout << "(countComponents...";
  }

  for (ulong i = 0; i < labels.rows(); i++) {
    for (ulong j = 0; j < labels.cols(); j++) {
      int k = labels.getPixel(i,j);
      if (k > 0) components[k].second += 1;
	}
  }
  sort(components.begin(), components.end(), greaterFrequency);
  
  int count = 0;
  for(ulong i = 0; i < components.size(); i++) {
    if (components[i].second > 0) {
      count++;
    } else {
      break;
    }
  }

  vector<int> key(components.size() + 1, -1);
  key[0] = 0;
  for(int i = 0; i < count; i++) {
    key[components[i].first] = i+1;
  }

  for (ulong i = 0; i < labels.rows(); i++) {
    for (ulong j = 0; j < labels.cols(); j++) {
      int k = labels.getPixel(i, j);
      int y = key[k];
      if (y == -1) {
        cerr << "countComponents(...): "
             << "Something's wrong in file " << __FILE__ << ", line " << __LINE__ << endl
             << "The sorted vector components appears to be corrupt." << endl;
        abort();
      }
      labels.setPixel(i, j, y);
    }
  }

  components.resize(count);
  if (verbosity > 1) {
     cout << "done)";
  }
  return count;
}

/* Used to sort the labels that are represented by pair.
 */
bool greaterFrequency(pair<int,int> lhs, pair<int,int> rhs) {
  return lhs.second > rhs.second;
}

// function labelComponents identifies the groups of pixels in the binary image input that are connected. Each pixel is then
// assigned a component label, which is returned through the second argument, labels. We assume that each interior pixel has
// exactly eight neighbors, while those on the boundary of the input image will have fewer.

int labelComponents(binaryImage &input, grayImage<GLint> &labels) {
  if (verbosity > 1) {
     cout << "labelComponents...";
  }
  int nextLabel = 1;
  labels.clear();
  labels.setRows(input.rows());
  labels.setCols(input.cols());
  labels.resize();
	
  for(ulong i = 0; i < input.rows(); i++)
    for (ulong j = 0; j < input.cols(); j++) {
      if (input.getPixel(i, j)) {
        for (ulong ii = i-1; ii <= i; ii++)
          for (ulong jj = j-1; jj < 2*(i - ii) + j; jj++) {
            int neighbor = labels.getPixel(ii, jj);						
            if (neighbor > 0) {
              int center = labels.getPixel(i, j);
              if (center == 0) {
                labels.setPixel(i, j, static_cast<GLint>(neighbor));
              } else if (neighbor != center) {
                mergeComponents(labels, i, j, max(neighbor, center), min(neighbor, center));
              }
            }
          }
        if (labels.getPixel(i, j) == 0) {                                // no upper-left neighboring pixels are set.
          labels.setPixel(i, j, static_cast<GLint>(nextLabel++));      // select a new component label.
        }
      }
	}
  if (verbosity > 1) {
     cout << "done\n";
  }
  return nextLabel;
}

// labelComponents identifies groups of pixels in the gray-level image input that are connected. Each pixel is then
// assigned a component label, which is returned through the second argument labels. The is an overloaded version of
// the previous version of this function.
int labelComponents(grayImage<GLubyte> &input, grayImage<GLint> &labels) {
  if (verbosity > 1) {
     cout << "labelComponents...";
  }
  int nextLabel = 1;
  labels.clear();
  labels.setRows(input.rows());
  labels.setCols(input.cols());
  labels.resize();
 
  for(ulong i = 0; i < input.rows(); i++)
    for (ulong j = 0; j < input.cols(); j++) {
      if (input.getPixel(i, j)) {
        for (ulong ii = i-1; ii <= i; ii++)
          for (ulong jj = j-1; jj < 2*(i - ii) + j; jj++) {
            int neighbor = labels.getPixel(ii, jj);						
            if (neighbor > 0) {
              int center = labels.getPixel(i, j);
              if (center == 0) {
                labels.setPixel(i, j, static_cast<GLint>(neighbor));
              } else if (neighbor != center) {
                mergeComponents(labels, i, j, max(neighbor, center), min(neighbor, center));
              }
            }
          }
        if (labels.getPixel(i, j) == 0) {                                // no upper-left neighboring pixels are set.
          labels.setPixel(i, j, static_cast<GLint>(nextLabel++));      // select a new component label.
        }
      }
	}
  if (verbosity > 1) {
     cout << "done\n";
  }
  return nextLabel;
}


// merge components assigns all pixel labels that have a value of from (argument 4), to the new value to (argument 5), in
// the indicated gray level image a, that contains the labels. Note that only pixels in the first r rows are affected, as
// are only those in columns 0 through c in the final row.
void mergeComponents(grayImage<GLint> &a, ulong r, ulong c, int from, int to) {
  assert(r < a.rows());
  assert(c < a.cols());
  for(ulong i = 0; i <= r; i++) {
    for (ulong j = 0; j < (i == r ? c : a.cols()); j++) {
      if (a.getPixel(i, j) == from) {
        a.setPixel(i, j, static_cast<GLint>(to));
      }
    }
  }
  return;
}


GLubyte pkRGB(int r, int g, int b) {
  return static_cast<GLubyte>(((((r & 7) << 3) | (g & 7)) << 3) | (b & 3));
}

void makeColorLabels(int nLabels) {
  vector<GLubyte> table(max(nLabels,8) +1, 0xff);
 
  table[0] = 0;             // black
  table[1] = pkRGB(7,0,0);  // red
  table[2] = pkRGB(5,4,0);  // orange
  table[3] = pkRGB(7,7,0);  // yellow
  table[4] = pkRGB(1,7,2);  // green
  table[5] = pkRGB(1,1,3);  // blue
  table[6] = pkRGB(0,1,2);  // magenta
  table[7] = pkRGB(5,2,3);  // magenta
  table[8] = pkRGB(3,0,3);  // magenta
 

  ulong rows = componentLabels.rows();
  ulong cols = componentLabels.cols();

  colorLabels.setRows(rows);
  colorLabels.setCols(cols);
  colorLabels.resize();
  for(ulong i = 0; i < rows; i++) {
    for(ulong j = 0; j < cols; j++) {
      ulong x = componentLabels.getPixel(i, j);
      colorLabels.setPixel(i, j, table[x]);
    }
  }
}



int collectComponents() {
  if (verbosity > 1) {
    cout << "(Computing component labels...";
  }
  closureOp();
  size_t nLabels = labelComponents(processedImage[1], componentLabels);

  if (verbosity > 1) {
    cout << "done.)" << endl;
    cout << "nLabels = " << nLabels << endl;
  }
 
  labelCount.clear();
  labelCount.resize(nLabels);
  for(size_t i = 0; i < nLabels; i++) {
    labelCount[i].first = i;
  }

  // Sort clusters by size, and count them.
  int clusterCount = countComponents(componentLabels, labelCount);

 if (verbosity > 1) {
   cout << endl
        << "Cluster count = " << clusterCount 
        << endl;

   for(int i = 0; i < min(30, clusterCount); i++) {
     cout << labelCount[i].second << " ";
   }
   cout << endl;
 }

  return clusterCount;
}

// computeComponentHull computes the two-dimensional convex hull of the indicated component.
void computeComponentHull(int label, ConvexHull2d<double> &h) {
  vector<ntuple<double,2> > points;

  for(ulong i = 0; i < componentLabels.rows(); i++) {
    for(ulong j = 0; j < componentLabels.cols(); j++) {
      if (componentLabels.getPixel(i, j) == label) {
        points.push_back(ntuple<double,2>(static_cast<double>(j), 
                                          static_cast<double>(i)));
      }
    }
  }

  h.compute(points);
}

// computeEigenvalues computes the eigenvalues and eigenvectors of the two-dimensional covariance matrix
// for the centered distribution of points associated with the component indexed by the value of label, the
// first parmater. Cluster is a struct defined in header project.h.
void computeEigenvalues(int label, Cluster &c) {
	double xm = 0;
	double ym = 0;
	double xym = 0;
	double x2m = 0;
	double y2m = 0;
	ulong count = 0;
	
	// For every pixel
	for(ulong i = 0; i < componentLabels.rows(); i++) {
	  for (ulong j = 0; j < componentLabels.cols(); j++) {
		if (componentLabels.getPixel(i, j) == label) {
		  xm += j;   //	Note that j indicates the column, i.e., the horizontal pixel offset.
		  ym += i;
		  x2m += j*j;
		  y2m += i*i;
		  xym += i*j;
		  count++;
		}
	  }
	}
	
	xm = xm/count;  	// mean x coordinate offset
	ym = ym/count;		// mean y coordiante offset
	x2m = x2m/count;    // mean x^2
	y2m = y2m/count;	// mean y^2
	xym = xym/count;	// mean x*y
	
	double sxx = x2m - xm*xm;						// <x^2> - <x>^2
	double sxy = xym - xm*ym;						// <x*y> - <x>*<y>
	double syy = y2m - ym*ym;						// <y^2> - <y>^2
	double splus = sxx + syy;						// <x^2> + <y^2> - (<x>^2 + <y>^2)
	double sdiff = syy - sxx;						// <y^2> - <x^2> - (<y>^2 - <x>^2)
	double sxy2 = sxy*sxy;
	double disc = sqrt(sdiff*sdiff + 4*sxy2);
	
	c.size = count;
	c.center[0] = xm;
	c.center[1] = ym;
	c.lam1 = (splus + disc)/2;  // the dominant eigenvalue
	c.lam2 = (splus - disc)/2;  // the subdominant eigenvalue

	// Solving the (homogeneous) linear system for the eigenvectors:
	//  x (sxx - lambda) = - y sxy
    //  x sxy = - (syy - lambda) y
    //
	// yields tan(theta) = y/x = - (sxx - lambda)/sxy, [which also equals - sxy/(syy - lambda)], whence,
	//
	// cos^2(theta) = 1/(1 + tan^2(theta)) = 2 sxy^2/((sxx-syy)^2 +/- (syy - sxx)sqrt((sxx-syy)^2 + 4 sxy^2) + 4 sxy^2),
	//
	// or,
	//
	// x = cos(theta) = sqrt(2 sxy^2 / (sdiff^2 + 4 sxy^2 +/- sdiff*sqrt(sdiff^2 + 4 sxy^2)))
	//   = sqrt(2 sxy2 / (sdiff^2 + 4 sxy2 +/- sdiff*disc))
    // 
	// and,
	//
	// y = sin(theta) = tan(theta) cos(theta) = x (sdiff +/- disc)/(2sxy)
	//
	double a2 = (2*sxy2/(sdiff*sdiff + sdiff*disc + 4*sxy2));
	double a  = sqrt(a2);  // x component of the dominant eigenvector
	double b2 = (2*sxy2/(sdiff*sdiff - sdiff*disc + 4*sxy2));
	double b  = sqrt(b2);  // x component of the subdominant eigenvector
		
	c.ev1[0] = a;
	c.ev1[1] = a*(sdiff + disc)/(2*sxy);
	c.ev2[0] = b;
	c.ev2[1] = b*(sdiff - disc)/(2*sxy);

	return;
}

