#ifndef __lsmImage_H__
#define __lsmImage_H__

//
//  lsmImage.h
//
//  Project: nuclearSurface.
//
//  Created by Robert R. Snapp on 2007-04-25.
//  Copyright (c)  Robert R. Snapp. All rights reserved.
//

#include "env.h"
#include "grayImage3d.h"
#include "lsm.h"
#include "lsmFileInfo.h"
//#include "imageComponentArray.h"
#include <iostream>

class LSM_Image {
public:
	explicit LSM_Image(LSM_File* lsm, int channelIndex = 0);
//	explicit LSM_Image(const char* filename, int channelIndex = 0);
	//void applyThreshold(int threshold);
	
	LSM_Image(const LSM_Image &im) : channel_(im.channel_), image_(im.image_) {
		info_ = new LSM_File_Info(*im.info_);
		lsmFile_ = new LSM_File(*im.lsmFile_);
	}
	
	LSM_Image& operator=(const LSM_Image &im) {
		// Guard against self-assignment
		if (this != &im) {
			if (info_) delete [] info_;
			info_ = new LSM_File_Info(*im.info_);
			if (lsmFile_) delete lsmFile_;
			lsmFile_ = new LSM_File(*im.lsmFile_);
			channel_ = im.channel_;
			image_   = im.image_;
		}
		return *this;
	}
	
	~LSM_Image() {
		delete info_;
		delete lsmFile_;
	}
	
	class bad_format{
		void debug_print() const {std::cerr << "bad file format error"; }
	};
	
	class bad_filename{
		void debug_print() const {std::cerr << "bad filename path error"; }
	};
	
//private:
	LSM_File_Info        *info_;
	LSM_File             *lsmFile_;
	grayImage3d<GLubyte> image_;
	int channel_;
	// ImageComponentArray *component_;
};


#endif // __lsmImage_H__
