#ifndef _LSM_FILE_INFO_H_
#define _LSM_FILE_INFO_H_

//
//  file: imageData.h
//  project: lsmModeler
// 
//  manages the data acquired from an LSM image.
//
//  Created by Robert Snapp on 2007-01-12.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "lsm.h"
#include "grayImage3d.h"

// Uses the interface functions of LSM_File to construct a concise inventory of the
// of the indicated lsm file. The number of images contained in the file is returned.
class LSM_File_Info {
private:
	struct directoryRecord {
		uint32 width;
		uint32 length;
		int depth;  // number of channels
		bool thumbnail;
	};
	uint32 info_size_;
	uint32 nOriginals_;
	uint32 nThumbnails_;
	directoryRecord *info_;

protected:
	static const bool debug = true;

public:
	static const int SEARCH_FAILED = -1;
	
	LSM_File_Info(LSM_File &lsm) : nOriginals_(0), nThumbnails_(0) {
		info_size_ = lsm.getImageCount();
		info_ = (directoryRecord*) new directoryRecord[info_size_];
		
		for(int i = 0; i < info_size_; i++) {
			info_[i].width     = lsm.getWidth(i);   // returns the (horizontal) width of the i-th file.
			info_[i].length    = lsm.getLength(i);  // returns the (vertical) length of the i-th file.
			info_[i].depth     = lsm.getSamplesPerPixel(i);
			if (info_[i].thumbnail = lsm.isLowResCopy(i)) {
				nThumbnails_++;
			} else {
				nOriginals_++;
			}		
		}
	}
	
	LSM_File_Info(const LSM_File_Info &a) : 
		info_size_(a.info_size_), nOriginals_(a.nOriginals_), nThumbnails_(a.nThumbnails_) {
			info_ = (directoryRecord*) new directoryRecord[info_size_];
			for (int i = 0; i < info_size_; i++) {
				info_[i].width = a.info_[i].width;
				info_[i].length = a.info_[i].length;
				info_[i].depth = a.info_[i].depth;
			}
	}
	
	LSM_File_Info& operator=(const LSM_File_Info &a) {
		//Guard against self-assignment
		if (this != &a) {
			info_size_ = a.info_size_;
			nOriginals_ = a.nOriginals_;
			nThumbnails_ = a.nThumbnails_;
			if (info_) delete [] info_;
			info_ = (directoryRecord*) new directoryRecord[info_size_];
			for (int i = 0; i < info_size_; i++) {
				info_[i].width  = a.info_[i].width;
				info_[i].length = a.info_[i].length;
				info_[i].depth  = a.info_[i].depth;
			}
		}
		return *this;
	}
	
	~LSM_File_Info() {
		if (info_) delete [] info_;
	}

	std::ostream& dump(std::ostream &os) const {
		for(int i = 0; i < info_size_; i++) {
			os << "Image " << std::setw(2) << i << " has " 
	       	   << std::setw(5)  << info_[i].width << "*" << std::setw(5) << info_[i].length << " pixels,"
	       	   << info_[i].depth << " channels.";
	    	if (info_[i].thumbnail) {
	 	  		os << " Is a thumbnail.";
			}
			return os << std::endl;
		}
	}

	uint32 numberOfImages() const {
			return info_size_;
	}
	
	uint32 numberOfOriginals() const {
			return nOriginals_;
	}
	uint32 numberOfThumbnails() const {
			return nThumbnails_;
	}
	
	uint32 width(int i) const {
		if (i < 0 | i >= info_size_) {
			throw std::out_of_range("ImageInfo::width(i) argument out of range.");
		}
		return info_[i].width;
	}
	uint32 length(int i) const {
		if (i < 0 | i >= info_size_) {
			throw std::out_of_range("ImageInfo::length(i) argument out of range.");
		}
		return info_[i].length;
	}
	uint32 depth(int i) const {
		if (i < 0 | i >= info_size_) {
			throw std::out_of_range("ImageInfo::depth(i) argument out of range.");
		}
		return info_[i].depth;
	}
	
	// findNextOriginalImage examines the indicated LSM_Directory, starting at position index,
	// for the next image that is not a thumbnail. The index is returned. If no such image
	// is found, the value SEARCH_FAILED is returned.
	int findNextOriginalImage(int index) const {
		for(int i = index; i < info_size_; i++) {
			if(! info_[i].thumbnail) return i;
		}
		return SEARCH_FAILED;
	}
	
	int findPreviousOriginalImage(int index) const {
		for(int i = index - 1; i >= 0; i--) {
			if(! info_[i].thumbnail) return i;
		}
		return SEARCH_FAILED;
	}

};


#endif /* _LSM_FILE_INFO_H_ */
