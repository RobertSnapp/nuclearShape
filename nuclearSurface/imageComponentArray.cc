//
//  imageComponentArray.cc
//  lsmModeler
//
//  Created by Robert Snapp on 2007-02-06.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//
#include "imageHistogram.h"
#include "nuclearSurface.h"
#include "imageComponentArray.h"
#include "grayImage3d.h"
#include "env.h"
#include "lsmVoxelImage.h"
#include "ntuple.h"
#include "convexhull.h"

#include <algorithm>
#include <stdexcept>
#include <set>

using namespace std;

// An STL pair<int,int> data structure is used to store a label in the first component,
// and its frequency in the second. Later we will need to sort the labels according to
// frequency; hence the following boolean function.
bool greaterFrequency(pair<int,int> lhs, pair<int,int> rhs) {
  	return lhs.second > rhs.second;
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

// computeComponents assigns a label to each voxel in the indicated three-dimensional binary
// image, such that voxels that are adjacent to one another (assuming a 26 point topology)
// are assigned a common label for that component, and the labels for each component are unique.
// This is a generalization of the operation described as "pixel labeling"
// Anil Jain's, <i>Fundamentals of Digital Image Processing</i>, Prentice-Hall, p. 409. 1989;
// and in Horn's <i>Robot Vision</i>. During the first pass, equivalence classes of labels are
// maintained in the vector parent, following an algorithm described in Knuth's <i>Art of Computer Programming</i>,
// Vol. 1, Second Edition, 1973, p. 353, attributed to M. J. Fischer and B. A. Galler.
void ImageComponentArray::computeComponents(const binaryImage3d &input) {
  if (debug) {
    cerr << "ImageComponentArray::computeComponents(...)" << endl;
  }

  // Initialization

  ulong nextLabel = 1;  // Start from 1 as 0 represents the background.

  ulong rows = input.rows();
  ulong cols = input.cols();
  ulong lays = input.layers();

  float dx = input.getDX();
  float dy = input.getDY();
  float dz = input.getDZ();
  
  d_labels.setRows(rows);
  d_labels.setCols(cols);
  d_labels.setLayers(lays);

  d_labels.setDX(dx);
  d_labels.setDY(dy);
  d_labels.setDZ(dz);
  if (verbosity > 1) {
    std::cerr << "Set voxel dimensions to (dx, dy, dz) = ("
              << d_labels.getDX() << ", "
              << d_labels.getDY() << ", "
              << d_labels.getDZ() << ")" << std::endl;
  } 

  d_labels.resize();
  
  // For maintaining equivalence classes of labels. We create a graph of 1000 disjoint nodes. 
  // (New nodes will be added as needed.) Each node is either a root node of a tree or is a successor 
  // in a tree. Using the vector<int> parent, we maintain this equivalence graph, as parent[n]
  // returns the parent of node n, i.e., node n is a successor of node parent[n]. Note that if
  // node n is a root node, then parent[n] equals nil.

  int const nil = -1;             // since vectors are zero-indexed.
  vector<int> parent(1000, nil); 
  

  if (verbosity > 1) {
    cerr << "First pass within computeComponents" << endl;
  }

  // First pass
  for(int l = 0; l < static_cast<int>(lays); l++) {
    if (verbosity > 1) {
      cerr << "L";
    }
    for(int i = 0; i < static_cast<int>(rows); i++) {
      for(int j = 0; j < static_cast<int>(cols); j++) {
        if(input.getVoxelClip(i, j, l)) {
          // For each active binary voxel, check if a positive label has been assigned to a previously
          // visited neighbor (assuming the 26-point topology). In general, this includes nine pixels
          // in the previous layer, three pixels in the previous row of the current layer, and the single
          // pixel in the previous column of the current row. However, the limits of the following for loops
          // take all special cases into account.
          for(int ll = max<int>(0, l-1); ll <= l; ll++)
            for(int ii = max<int>(0, i-1); ii < min<int>(rows, i + l - ll + 1); ii++)
              for(int jj = max<int>(0, j-1); jj < min<int>(cols, ((ll < l || ii < i) ? j + 2 : j)); jj++) {
                // read the label of the neighboring voxel
                int neighbor = d_labels.getVoxel(ii, jj, ll);
                if (neighbor > 0) {
                  // A positive value indicates that the neighboring pixel is active, and has been assigned
                  // to the corresponding component. 
                  int center = d_labels.getVoxel(i, j, l);
                  if (center == 0) {
                    // A zero value indicates that the current voxel has not yet been assigned to a component.
                    // Thus, we assign it to the same component as its neighbor.
                    d_labels.setVoxel(i, j, l, static_cast<ulong>(neighbor));
                  } else if (neighbor != center) {
                    // A different positive value indicates that this voxel unites two components that were
                    // previously assumed to be disjoint. We merge the two components by placing their labels
                    // in a common equivalence class. (See Knuth, vol. 1).
                    while(parent[center] != nil) {
                      center = parent[center];
                    }

                    while(parent[neighbor] != nil) {
                      neighbor = parent[neighbor];
                    }

                    if(center != neighbor) {
                      parent[neighbor] = center;
                    }
                  }
                }
              }

          if (d_labels.getVoxel(i, j, l) == 0) {                            
            // No upper-left neighboring pixels are set. Select a new component label.
            d_labels.setVoxel(i, j, l, static_cast<GLint>(nextLabel++));

            if (nextLabel >= parent.size()) {
              // If necessary, enlarge the equivalence graph, parent.
              parent.resize(2*nextLabel, nil);
            }   
          } 
        }   
      } 
    }   
  }	    

  if (verbosity > 1) {
	cerr << endl
		 << nextLabel - 1 << " labels generated."
		 << endl;
  }

  // Consolidate the equivalence graph: replace the value of parent[n] with the index of the
  // root node of the tree that contains node n.
  for(ulong i = 1; i < nextLabel; ++i) {
    if(parent[i] != nil) {
      // If i is not a root node:
      int rootNode;
      while((rootNode = parent[parent[i]]) != nil) {
        parent[i] = rootNode;
      }
    }
  }

  if (verbosity > 1) {
    cerr << "\nSecond pass." << endl;
  }

  // Second pass: Replace each label with that of its root node.
  for(ulong l = 0; l < lays; l++) {
    if (verbosity > 1) {
      cerr << "L";
    }
    for(ulong i = 0; i < rows; i++)
      for(ulong j = 0; j < cols; j++) {
        ulong k = d_labels.getVoxel(i, j, l);
        int v = parent[k];
        if (v != nil) {
          d_labels.setVoxel(i, j, l, static_cast<ulong>(v));
        }
      }
  }

   
  if (verbosity > 1) {
    cerr << "\nRank the positively labeled components." << endl;
  }
  // Rank the positively labeled components by frequency. 
  d_labelCount.resize(nextLabel);
  for(ulong k = 0; k < nextLabel; k++) {
    d_labelCount[k].first = k;
  }
	
  for (ulong l = 0; l < lays; l++)
    for (ulong i = 0; i < rows; i++)
      for (ulong j = 0; j < cols; j++) {
        ulong k = d_labels.getVoxel(i, j, l);
        if (k > 0) d_labelCount[k].second++;  // Ignore 0, which represents the background
	}
	
  std::sort(d_labelCount.begin(), d_labelCount.end(), greaterFrequency);

  if (verbosity > 1) {
	ulong const maxLabelsShown = 25;
    for(ulong i = 0; i < min(maxLabelsShown, nextLabel); i++) {
      cerr << d_labelCount[i].second << ", ";
    }
	cerr << endl;
  }

  //////
  // Compute the number of distinct labels in use.
 int count = 0;
 for(ulong i = 0; i < d_labelCount.size(); i++) {
   if (d_labelCount[i].second > 0) {
	 count++;
   } else {
	 break;
   }
 }

 // Convert the labels so that they assume the sequence 1, 2, 3, ..., with 1 identifying
 // the largest componet, 2 the second largest, etc., and 0 denotes the background. The
 // element key[i] will provide the value of the new label for what was formerly known as i.

 vector<int> key(d_labelCount.size() + 1, -1);  // -1 denotes "undefinded."
 key[0] = 0; // background label.
 for(int i = 0; i < count; i++) {
   key[d_labelCount[i].first] = i+1;
   d_labelCount[i].first = i+1;
 }

 // Assign a new label to each voxel, using the new scheme. Note that only the first count
 // elements of key have a value other than -1.
 for (ulong l = 0; l < d_labels.layers(); l++) {
   for (ulong i = 0; i < d_labels.rows(); i++) {
	 for (ulong j = 0; j < d_labels.cols(); j++) {
	   int k = d_labels.getVoxel(i, j, l);
	   int y = key[k];
	   if (y == -1) {
		 cerr << "computeComponents(...): "
			  << "Something's wrong in file " << __FILE__ << ", line " << __LINE__ << endl
			  << "The sorted vector components appears to be corrupt." << endl;
		 abort();
	   }
	   d_labels.setVoxel(i, j, l, y);
    }
  }

  d_labelCount.resize(count);
 }
}

#ifdef COMMENT
void ImageComponentArray::computeComponents(const grayImage3d<GLubyte> &input) {
  if (debug) {
    cerr << "ImageComponentArray::computeComponents(...)" << endl;
  }

  // Initialization

  ulong nextLabel = 1;  // Start from 1 as 0 represents the background.

  ulong rows = input.rows();
  ulong cols = input.cols();
  ulong lays = input.layers();

  float dx = input.getDX();
  float dy = input.getDY();
  float dz = input.getDZ();
  
  d_labels.setRows(rows);
  d_labels.setCols(cols);
  d_labels.setLayers(lays);

  d_labels.setDX(dx);
  d_labels.setDY(dy);
  d_labels.setDZ(dz);
  if (verbosity > 1) {
    std::cerr << "Set voxel dimensions to (dx, dy, dz) = ("
              << d_labels.getDX() << ", "
              << d_labels.getDY() << ", "
              << d_labels.getDZ() << ")" << std::endl;
  } 

  d_labels.resize();

  // For maintaining equivalence classes of labels. We create a graph of 1000 disjoint nodes. 
  // (New nodes will be added as needed.) Each node is either a root node of a tree or is a successor 
  // in a tree. Using the vector<int> parent, we maintain this equivalence graph, as parent[n]
  // returns the parent of node n, i.e., node n is a successor of node parent[n]. Note that if
  // node n is a root node, then parent[n] equals nil.

  int const nil = -1;             // since vectors are zero-indexed.
  vector<int> parent(1000, nil); 
  

  if (verbosity > 1) {
    cerr << "First pass within computeComponents" << endl;
  }

  // First pass
  for(int l = 0; l < static_cast<int>(lays); l++) {
    if (verbosity > 1) {
      cerr << "L";
    }
    for(int i = 0; i < static_cast<int>(rows); i++) {
      for(int j = 0; j < static_cast<int>(cols); j++) {
        if(input.getVoxelClip(i, j, l)) {
          // For each active binary voxel, check if a positive label has been assigned to a previously
          // visited neighbor (assuming the 26-point topology). In general, this includes nine pixels
          // in the previous layer, three pixels in the previous row of the current layer, and the single
          // pixel in the previous column of the current row. However, the limits of the following for loops
          // take all special cases into account.
          for(int ll = max<int>(0, l-1); ll <= l; ll++)
            for(int ii = max<int>(0, i-1); ii < min<int>(rows, i + l - ll + 1); ii++)
              for(int jj = max<int>(0, j-1); jj < min<int>(cols, ((ll < l || ii < i) ? j + 2 : j)); jj++) {
                // read the label of the neighboring voxel
                int neighbor = d_labels.getVoxel(ii, jj, ll);
                if (neighbor > 0) {
                  // A positive value indicates that the neighboring pixel is active, and has been assigned
                  // to the corresponding component. 
                  int center = d_labels.getVoxel(i, j, l);
                  if (center == 0) {
                    // A zero value indicates that the current voxel has not yet been assigned to a component.
                    // Thus, we assign it to the same component as its neighbor.
                    d_labels.setVoxel(i, j, l, static_cast<ulong>(neighbor));
                  } else if (neighbor != center) {
                    // A different positive value indicates that this voxel unites two components that were
                    // previously assumed to be disjoint. We merge the two components by placing their labels
                    // in a common equivalence class. (See Knuth, vol. 1).
                    while(parent[center] != nil) {
                      center = parent[center];
                    }

                    while(parent[neighbor] != nil) {
                      neighbor = parent[neighbor];
                    }

                    if(center != neighbor) {
                      parent[neighbor] = center;
                    }
                  }
                }
              }

          if (d_labels.getVoxel(i, j, l) == 0) {                            
            // No upper-left neighboring pixels are set. Select a new component label.
            d_labels.setVoxel(i, j, l, static_cast<GLint>(nextLabel++));

            if (nextLabel >= parent.size()) {
              // If necessary, enlarge the equivalence graph, parent.
              parent.resize(2*nextLabel, nil);
            }   
          } 
        }   
      } 
    }   
  }	    

  if (verbosity > 1) {
	cerr << endl
		 << nextLabel << " labels generated."
		 << endl;
  }

  // Consolidate the equivalence graph: replace the value of parent[n] with the index of the
  // root node of the tree that contains node n.
  for(ulong i = 1; i < min(nextLabel, parent.size()); ++i) {
    if(parent[i] != nil) {
      // If i is not a root node:
      int rootNode;
      while((rootNode = parent[parent[i]]) != nil) {
        parent[i] = rootNode;
      }
    }
  }

  if (verbosity > 1) {
    cerr << "\nSecond pass." << endl;
  }

  // Second pass: Replace each label with that of its root node.
  for(ulong l = 0; l < lays; l++) {
    if (verbosity > 1) {
      cerr << "L";
    }
    for(ulong i = 0; i < rows; i++)
      for(ulong j = 0; j < cols; j++) {
        ulong k = d_labels.getVoxel(i, j, l);
        int v = parent[k];
        if (v != nil) {
          d_labels.setVoxel(i, j, l, static_cast<ulong>(v));
        }
      }
  }

   
  if (verbosity > 1) {
    cerr << "\nRank the positively labeled components." << endl;
  }
  // Rank the positively labeled components by frequency. 
  d_labelCount.resize(nextLabel);
  for(ulong k = 0; k < nextLabel; k++) {
    d_labelCount[k].first = k;
  }
	
  for (ulong l = 0; l < lays; l++)
    for (ulong i = 0; i < rows; i++)
      for (ulong j = 0; j < cols; j++) {
        ulong k = d_labels.getVoxel(i, j, l);
        if (k > 0) d_labelCount[k].second += 1;  // Ignore 0, which represents the background
	}
	
  std::sort(d_labelCount.begin(), d_labelCount.end(), greaterFrequency);

  if (debug) {
	ulong const maxLabelsShown = 25;
    for(ulong i = 0; i < min(maxLabelsShown, nextLabel); i++) {
      cerr << d_labelCount[i].second << ", ";
    }
  }
	
  if (debug) cerr << "done!" << endl;
}

#endif

#ifdef COMMENT
// The following constructor, ImageComponentArray(...),
// creates an image component array directly from the indicated channel of an lsm image.
// All voxels with intensities greater than or equal to the specified threshold will
// be classified as being active. The domain (or field of view) of the voxel volume
// can be restricted by specifing an appropriate subvolume.

ImageComponentArray::ImageComponentArray(const lsmVoxelImage &input, 
										ulong channel,
                                        Subvolume &subvol,
										GLubyte threshold) {
	if (debug) {
		cerr << "ImageComponentArray(...):" << endl;
	}
	
	ulong nextLabel = 1;
	uint32 rows   = subvol.height; // input.rows();
	uint32 cols   = subvol.width;  // input.cols();
	uint32 layers = subvol.depth;  // input.layers();
	
	float dx = input.getDX();
	float dy = input.getDY();
	float dz = input.getDZ();
	
	if (threshold == 0) {
		// select a threshold using the following adaptive algorithm
		ImageHistogram h(1);
		h.processImage(input, channel);
		threshold = static_cast<ulong>(h.getMean() + 1.5*h.getSDev());	
	}
	
	if(debug){
		cerr << "rows = " << rows << endl
		     << "cols = " << cols << endl
		     << "layers = " << layers << endl
		     << "dx = " << dx << " meters" << endl
			 << "dy = " << dy << " meters" << endl
			 << "dz = " << dz << " meters" << endl
	    	 <<"actual threshold = " << threshold << endl;
	}
	
	d_labels.setDX(dx);
	d_labels.setDY(dy);
	d_labels.setDZ(dz);
	d_labels.setRows(rows);
	d_labels.setCols(cols);
	d_labels.setLayers(layers);
	d_labels.resize();

    int const root = -1;

    vector<int> parent(1000, root);

    for(int l = 0; l < static_cast<int>(layers); l++) {
      if (debug) cerr << "!"; 						// poor man's progress bar.
      int lImage = l + subvol.layer;
      for(int i = 0; i < static_cast<int>(rows); i++) {
        int iImage = i + subvol.row;
        if (debug && i % 128 == 0) cerr << ".";		// poor man's progress bar.
        for(int j = 0; j < static_cast<int>(cols); j++) {
          int jImage = j + subvol.col;
          if (input.getVoxelClip(iImage, jImage, lImage, channel) >= threshold) {
            //	Assign a label to the current voxel that matches the label of any contiguous
            //  neighbors, assuming a 26-point topology, of which at most 9 + 3 + 1 = 13
            //  adjacent neighbors have been examined.
				
            for(int ll = max(0, l-1); ll <= l; ll++) 
              for (int ii = max<int>(0, i-1); ii <  min<int>(rows, i + l - ll + 1); ii++)
                for (int jj = max<int>(0, j-1); jj < min<int>(cols, ((ll < l || ii < i) ? j + 2 : j)); jj++) {
                  int neighbor = d_labels.getVoxel(ii, jj, ll);
                  if (neighbor > 0) {
                    int center = d_labels.getVoxel(i, j, l);
                    if (center == 0) {
                      d_labels.setVoxel(i, j, l, static_cast<ulong>(neighbor));
                    } else if (neighbor != center) {
                      int maxValue = max<int>(neighbor, center);
                      
                      if (maxValue >= static_cast<int>(parent.size())) {
                        parent.resize(2*maxValue, root);
                      }

                      while(parent[center] != root) {
                        center = parent[center];
                      }

                      while(parent[neighbor] != root) {
                        neighbor = parent[neighbor];
                      }

                      if(center != neighbor) {
                        parent[neighbor] = center;
                      }
                    }
                  }
                }
					
            if (d_labels.getVoxel(i, j, l) == 0) {                            
              // no upper-left neighboring pixels are set. Select a new component label.
              d_labels.setVoxel(i, j, l, static_cast<GLint>(nextLabel++));      
            }
          }
        }
      }
	}	

	if (debug) cerr << nextLabel << " labels generated, ";


    // Make the equivalence graph shallow.
    for(ulong i = 1; i < min(nextLabel, parent.size()); ++i) {
      if(parent[i] != root) {
        int gramps;
        while((gramps = parent[parent[i]]) != root) {
          parent[i] = gramps;
        }
      }
    }

    for(ulong l = 0; l < layers; l++)
      for(ulong i = 0; i < rows; i++)
        for(ulong j = 0; j < cols; j++) {
          ulong k = d_labels.getVoxel(i, j, l);
          int v = parent[k];
          if (v != root) {
            d_labels.setVoxel(i, j, l, v);
          }
        }

	// Order the positively labeled components, according to their prominence. 
	d_labelCount.resize(nextLabel);
	for(ulong k = 0; k < nextLabel; k++) {
		d_labelCount[k].first = k;
	}
	
	for (ulong l = 0; l < layers; l++)
		for (ulong i = 0; i < rows; i++)
			for (ulong j = 0; j < cols; j++) {
				ulong k = d_labels.getVoxel(i, j, l);
				if (k > 0) d_labelCount[k].second += 1;
	}
	
	std::sort(d_labelCount.begin(), d_labelCount.end(), greaterFrequency);

	if (debug) {
	  ulong const maxLabelsShown = 25;
	  for(ulong i = 0; i < min(maxLabelsShown, nextLabel); i++) {
		cerr << d_labelCount[i].second << ", ";
	  }
	}
	
	if (debug) cerr << "done!" << endl;
}
#endif




// This version of pushComponentVoxels pushes only the coordinates of the voxel centers onto 
// the vector container.
void ImageComponentArray::pushComponentVoxels(int rank, std::vector<ntuple<int, 3> > &points) {
  if (rank < 0 || rank > static_cast<int>(d_labelCount.size()))
	throw std::out_of_range("pushComponentVoxels: bad value for rank.");
	
  ulong label = d_labelCount[rank].first;
  for(ulong l = 0; l < d_labels.layers(); ++l)
	for(ulong i = 0; i < d_labels.rows(); ++i)
	  for(ulong j = 0; j < d_labels.cols(); ++j) {
		if (d_labels.getVoxel(i, j, l) == label) {
		  points.push_back(ntuple<int, 3>(static_cast<int>(i), 
										  static_cast<int>(j), 
										  static_cast<int>(l)));
		}
	  }
  return;
}


#ifdef COMMENT
// This version of pushComponentVoxels inserts the coordinates of the eight vertices that
// delimit each active voxel onto the vector points.
void ImageComponentArray::pushComponentVoxels(int rank, std::vector<ntuple<int, 3> > &points) {
  if (rank < 0 || rank > static_cast<int>(d_labelCount.size()))
	throw std::out_of_range("pushComponentVoxels: bad value for rank.");
	
  ulong label = d_labelCount[rank].first;
  for(ulong l = 0; l < d_labels.layers(); ++l)
	for(ulong i = 0; i < d_labels.rows(); ++i)
	  for(ulong j = 0; j < d_labels.cols(); ++j) {
		if (d_labels.getVoxel(i, j, l) == label) {
		  for(ulong ll = l; ll <= l + 1; ll++) {
			for(ulong ii = i; ii <= i + 1; ii++) {
			  for(ulong jj = j; jj <= j + 1; jj++) {
				points.push_back(ntuple<int, 3>(static_cast<int>(ii), 
												static_cast<int>(jj), 
												static_cast<int>(ll)));
			  }
			}
		  }
		}
	  }

  return;
}
#endif



void ImageComponentArray::pushComponentShell(ulong rank, std::vector<ntuple<int, 3> > &points) {
  if (rank > d_labelCount.size()) {
     throw std::out_of_range("pushComponentShell: bad value for rank.");
  }
  ulong layers = d_labels.layers();
  ulong rows   = d_labels.rows();
  ulong cols   = d_labels.cols();
  ulong label  = d_labelCount[rank].first;

  for(int l = 0; l < static_cast<int>(layers); ++l)
	for(int i = 0; i < static_cast<int>(rows); ++i)
	  for(int j = 0; j < static_cast<int>(cols); ++j) {
              if (d_labels.getVoxel(i, j, l) == label) {
                // Examine the 26 possible neighbors of voxel (i, j, l)
                bool isShellVoxel = false;
                for (int ll = max(0, l - 1); ll < min(static_cast<int>(layers), l + 1); ll++)
                  for (int ii = max(0, i - 1); ii < min(static_cast<int>(rows), i + 1); ii++)
                    for (int jj = max(0, j - 1); jj < min(static_cast<int>(cols), j + 1); jj++) {
                      if (d_labels.getVoxel(ii, jj, ll) != label) {
                        isShellVoxel = true;
                        break;
                      }
                    }
                if (isShellVoxel) {
                  for(int ll = l; ll <= l + 1; ll++) {
                    for(int ii = i; ii <= i + 1; ii++) {
                      for(int jj = j; jj <= j + 1; jj++) {
						points.push_back(ntuple<int, 3>(static_cast<int>(ii), 
														static_cast<int>(jj), 
														static_cast<int>(ll)));
                      }
                    }
                  }
				}
              }
	  }
  return;
}


void ImageComponentArray::renderVoxel(int i, int j, int l) {
	float dx = d_labels.getDX();
	float dy = d_labels.getDY();
	float dz = d_labels.getDZ();
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
	ulong targetLabel = d_labelCount[component].first;
	ulong voxelCount = 0;
	for(ulong l = 0; l < d_labels.layers(); l++)
		for(ulong i = 0; i < d_labels.rows(); i++)
			for(ulong j = 0; j < d_labels.cols(); j++) {
				if (targetLabel == d_labels.getVoxel(i, j, l)) {
					renderVoxel(i, j, l);
					voxelCount++;
				}
			}
			
	if (debug) cerr << voxelCount << " voxels rendered ... done!" << endl;
	return;
}

void ImageComponentArray::displayLayer(ulong l) {
  d_colors.rgbPixels(l);
}
  
double ImageComponentArray::volume(ulong component) const {
  double totalVol = 0;
  double voxelVol = d_labels.getDX()*d_labels.getDY()*d_labels.getDZ();
  ulong layers = d_labels.layers();
  ulong rows   = d_labels.rows();
  ulong cols   = d_labels.cols();
	
  ulong targetLabel = d_labelCount[component].first;
  for(int l = 0; l < static_cast<int>(layers); l++)
	for(int i = 0; i < static_cast<int>(rows); i++)
	  for(int j = 0; j < static_cast<int>(cols); j++) {
		if (d_labels.getVoxelClip(i, j, l) == targetLabel) {
		  int totalNeighbors = 0;
		  int likeNeighbors = 0;
		  for(int ll = max<int>(0, l-1); ll < min<int>(static_cast<int>(layers), l+2); ll++)
			for(int ii = max<int>(0, i-1); ii < min<int>(static_cast<int>(rows), i+2); ii++)
			  for(int jj = max<int>(0, j-1); jj < min<int>(static_cast<int>(cols), j+2); jj++) {
				totalNeighbors++;
				if (d_labels.getVoxelClip(ii, jj, ll) == targetLabel)
				  likeNeighbors++;
			  }
		  totalVol += (1 + static_cast<double>(likeNeighbors))/(1 + static_cast<double>(totalNeighbors));
		}
	  }
  return totalVol*voxelVol;			
}



ulong ImageComponentArray::boundaryVoxelCount(int component) const {
	ulong boundaryVoxels = 0;
	
	ulong layers = d_labels.layers();
	ulong rows   = d_labels.rows();
	ulong cols   = d_labels.cols();
	
	ulong targetLabel = d_labelCount[component].first;
	for(int l = 0; l < static_cast<int>(layers); l++)
	  for(int i = 0; i < static_cast<int>(rows); i++)
		for(int j = 0; j < static_cast<int>(cols); j++) {
		  if (d_labels.getVoxelClip(i, j, l) == targetLabel) {
			bool isBoundaryVoxel = false;
			for(int ll = max<int>(0, l-1); ll < min<int>(static_cast<int>(layers), l+2); ll++)
			  for(int ii = max<int>(0, i-1); ii < min<int>(static_cast<int>(rows), i+2); ii++)
				for(int jj = max<int>(0, j-1); jj < min<int>(static_cast<int>(cols), j+2); jj++) {
				  if (d_labels.getVoxelClip(ii, jj, ll) != targetLabel)
					isBoundaryVoxel = true;
				}
			if (isBoundaryVoxel) boundaryVoxels++;
		  }
		}
	return boundaryVoxels;			
}

void ImageComponentArray::computeConvexHull(ulong i, ConvexHull<int> &c, HullType mode) {
  vector<ntuple<int,3> > v;
  
  if (verbosity > 1) {
	cerr << "computeConvexHull: Component " << i 
		 << " uses integer label " << d_labelCount[i].first 
		 << ", and should have " << d_labelCount[i].second << " voxels." << endl;
  }

  switch (mode) {
  case VOXEL_CENTERS:
    pushComponentVoxels(i, v);
    break;
  case VOXEL_CORNERS:
    pushComponentShell(i, v);
    break;
  }

  if (verbosity > 2) {
	cout << "Creating a convex hull from " << v.size() << " points." << endl;
  }
  c.compute(v);
  return;
}


GLubyte pkRGB(int r, int g, int b) {
  return static_cast<GLubyte>(((((r & 7) << 3) | (g & 7)) << 3) | (b & 3));
}

void ImageComponentArray::makeColorLabels(size_t nLabels) {
  vector<GLubyte> table(max(nLabels, static_cast<size_t>(8)) + 1, 0xff);
 
  table[0] = 0;             // black
  table[1] = pkRGB(7,0,0);  // red
  table[2] = pkRGB(5,4,0);  // orange
  table[3] = pkRGB(7,7,0);  // yellow
  table[4] = pkRGB(1,7,2);  // green
  table[5] = pkRGB(1,1,3);  // blue
  table[6] = pkRGB(0,1,2);  // magenta
  table[7] = pkRGB(5,2,3);  // magenta
  table[8] = pkRGB(3,0,3);  // magenta
 

  ulong rows   = d_labels.rows();
  ulong cols   = d_labels.cols();
  ulong layers = d_labels.layers();

  d_colors.setRows(rows);
  d_colors.setCols(cols);
  d_colors.setLayers(layers);
  d_colors.resize();
  for(ulong l = 0; l < layers; l++) {
	for(ulong i = 0; i < rows; i++) {
	  for(ulong j = 0; j < cols; j++) {
		ulong x = d_labels.getVoxel(i, j, l);
		d_colors.setVoxel(i, j, l, table[x]);
	  }
    }
  }
}

