
/******************************************************************************
 *
 *	File:			main.cc
 *
 *	Function:		main function for testing liblsm
 *
 *	Author:			Robert R. Snapp
 *
 *	Copyright:		Copyright (C) 2004 Robert R. Snapp
 *					All Rights Reserved.
 *
 *	Source:			Original.
 *
 *	Notes:			
 *
 *	Change History:
 *			2006_12_25	Started source.
 *	
 ******************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include "lsm.h"
#include "env.h"
#include "lsmVoxelImage.h"


#if (DEBUG > 0)
#define debug(exp) exp
#else
#define debug(exp)
#endif


using namespace std;

/* eraseDirectories is a string utility that removes the leading subirectories 
 * from the string referenced by path. Thus, the following code:
 *
 * string path("/home/users/smith/data.dat");
 * eraseDirectories(path);
 * cout << path << endl;
 *
 * should send the string "data.dat" to standard output.
 */

string& eraseDirectories(string &path) {
 	size_t last_slash = path.find_last_of("/");  // Remove leading directories
   	if (last_slash != string::npos)
   		(void) path.erase(0, last_slash + 1);
   	return path;
}

int main (int argc, char * const argv[]) {

   	debug(cerr << "argc = " << argc << endl); 	
   	if (argc < 2) {
   		string command = argv[0];                       // Identify the command name
   		(void) eraseDirectories(command);  // Remove leading directories (if present).
   		 
   		cerr << "Usage: " << command << " <filename>" << endl;
   		return EXIT_FAILURE;
   	}
   
   	LSM_File lsm(argv[1]);
	int imageCount = lsm.getImageCount();
	
	cout << "File " << argv[1] << " contains " << imageCount << " images." << endl;

	for(int i = 0; i < imageCount; i++) {
		int width = lsm.getWidth(i);
		int length = lsm.getLength(i);
		int bands  = lsm.getSamplesPerPixel(i);
		bool thumbnail = lsm.isThumbnail(i);
		cout << "Image " << setw(2) << i << " has " 
		     << setw(5)  << width << "*" << setw(5) << length << " pixels,"
		     << bands << " channels.";
		if (thumbnail) {
			cout << " Is a thumbnail.";
		}
		cout << endl;
	}

	//#ifdef COMMENT	
	cout << "\n\nTesting iterator:\n";
	
   	int j = 0;
	LSM_File::iterator iter;
	LSM_File::iterator lsmEnd = lsm.end();
	for(iter = lsm.begin(); iter != lsmEnd; ++iter) {
	  cout << "offset(iter) = " << distance(lsm.begin(),iter) << ", ";
	  cout << "sizeof(*iter) = " << sizeof(*iter) << endl;
	  // iter->dump(cout);

	  int width = iter->getWidth();
	  int length = iter->getLength();
	  int bands  = iter->getSamplesPerPixel();
	  bool thumbnail = iter->isThumbnail();
	  cout << "Image " << setw(2) << j++ << " has " 
		   << setw(5)  << width << "*" << setw(5) << length << " pixels,"
		   << bands << " channels";
	  if (thumbnail) {
		cout << ", (thumbnail).";
	  } else {
		cout << ".";
	  }
	  cout << endl;

	}
	//#endif

	cout << "\n\nTesting operator[]:\n";
	
	
	for(int i = 0; i < imageCount; ++i) {
		int width      = lsm[i].getWidth();
		int length     = lsm[i].getLength();
		int bands      = lsm[i].getSamplesPerPixel();
		bool thumbnail = lsm[i].isThumbnail();
		cout << "Image " << setw(2) << i << " has " 
		     << setw(5)  << width << "*" << setw(5) << length << " pixels,"
		     << bands << " channels.";
		if (thumbnail) {
			cout << " Is a thumbnail.";
		}
		cout << endl;
	}
	
	lsm.close();
	
	cout << "\n\nTesting lsmVoxelImage.h:\n";
	
	cout << "Instantiating lsmVoxelImage with channel = 0: lsmVoxelImage image1(argv[1], 0) with argv[1] = " 
	     << argv[1] << "...";
	lsmVoxelImage image1(argv[1], 0);
	cout << "done." << endl;
	
	cout << "rows   = " << image1.rows() << endl
	     << "cols   = " << image1.cols() << endl
		 << "bands  = " << image1.bands() << endl
		 << "layers = " << image1.layers() << endl;
		
	
	image1.clear(); // Clear the image, close the input stream.

	cout << "Instantiating lsmVoxelImage with all channels: lsmVoxelImage image2(argv[1]) with argv[1] = " 
	     << argv[1] << "...";
	lsmVoxelImage image2(argv[1]);
	cout << "done." << endl;

	cout << "rows   = " << image2.rows() << endl
	     << "cols   = " << image2.cols() << endl
		 << "bands  = " << image2.bands() << endl
		 << "layers = " << image2.layers() << endl;
	
	// return 0;
}
