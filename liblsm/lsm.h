/*
 *  lsm.h
 *  Project liblsm, a C++ library for reading and parsing data files produced
 *  by Zeiss laser scanning confocal microscopes. These files typically have
 *  an extension ".lsm" and contain a sequence of multichannel images with
 *  low resolution "thumbnails."
 *
 *  Usage: To read and parse the lsm data in filename.lsm, invoke
 *
 * 		#include "lsm.h"
 * 		...
 * 		LSM_File dataFile("filename.lsm");
 *		std::cout << "Found " << dataFile.getImageCount() << " images.\n";
 *
 *  Data for individual images in the series can be accessed by three different
 *  interfaces:
 *  
 *  by using a set of member functions, e.g. dataFile.getLength(int i) returns
 *  the number of rows in the i-th image in the series;
 *
 *  by using an [] operator, e.g., dataFile[i].getLength(), returns the number
 *  of rows of the same image; or
 *
 *  by using an iterator:
 *
 *     for(LSM_File::iterator iter = dataFile.begin(); 
 *         iter != dataFile.end(); ++iter) {
 *       cout << iter->getLength() << "x" << iter->getWidth() << endl;
 *     }
 *
 *  Created by Robert Snapp on 1/4/07.
 *  Copyright 2007 Robert R. Snapp. All rights reserved.
 *
 */

#ifndef _LSM_H_
#define _LSM_H_

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tiff.h>  // Defines uint8, uint16, uint32, and TIFFTAG_s.

typedef signed long sint32;
typedef double float64;

typedef struct LSM_Color_t {
	uint8 r;
	uint8 g;
	uint8 b;
} LSM_Color;

#if (DEBUG > 5)
#define lsm_debug(exp) exp
#else
#define lsm_debug(exp)
#endif

class LSM_File {
public:
	static const uint16 LSM_MAGIC_NUMBER = 0x2a;

	// Constructor: Given a path to an LSM file, the constructor parses the
	// header and parses each LSM subdirectory.
	LSM_File(const char* path);
	
	// Destructor: closes the file stream.
	~LSM_File();
	
	// Assume for now the default copy constructor, etc.
	
	// Get the location of the first directory in the current LMS file.
	uint32 getFirstOffset() {
		return d_firstDirectoryOffset;
	}
	
	void close() { if (d_in->is_open()) d_in->close(); }
	
	// Get the location of the cz_lmsinfo structure in the current LMS file.
	void parseCZInfo(uint32 offset);
	
	void parseColorTable(uint32 offset);
	
	void parseScanInfo(uint32 offset);
	
	// getImageCount returns the number of images found in the current file.
	uint32 getImageCount(void) const;

	// isThumbnail returns true if the indicated image is a low 
	// resolution copy (or thumbnail) of another image within the file. 
	// A value of false indicates that this image is full resolution.
	bool isThumbnail(uint32 image) const;
	
	// isPlanar returns true if the indicated image is in "planar format",
	// in which each component is stored as a separate plane. A value of 
	// false is returned if the components are stored in "chunky format",
	// in which the components are interleaved: RGBRGBRGB ...
	// Note that if there is only one component, there is no distinction
	// between the two. In this case, isPlanar also returns true.
	bool isPlanar(uint32 image) const;
	
	// isRGB returns true if the indicated image is stored in RGB format.
	bool isRGB(uint32 image) const;

	// getStripByteCount returns the number of bytes in the indicated image
	// channel for the indicated image.
	uint32 getStripByteCount(uint32 image, uint32 channel) const;

	// getStripOffset returns the location of the image data for then indicated
	// channel in the indicated image.
	uint32 getStripOffset(uint32 image, uint32 channel) const;

	// getSamplesPerPixel(void) returns the number of channels in the
	// indicated image.
	uint16 getSamplesPerPixel(uint32 image) const;

	// getBitsPerSample returns the number of bits for the indicated channel
	// in the indicated image. This value is typically 8 for each channel of
	// 24-bit RGB image, for which channel should equal 0, 1, or 2.
	uint16 getBitsPerSample(uint32 image, uint32 channel) const;
		
	// getLength returns the (vertical) length of the indicated image in pixels.
	uint32 getHeight(uint32 image) const;
	uint32 getLength(uint32 image) const {return getHeight(image);}
	
	// getWidth returns the (horizontal) width of the indicated image in pixels.
	uint32 getWidth(uint32 image) const;
	
	// getPixels copies a block of image data from the current LSM file, to a
	// block of memory indicated by the third argument. The number of bytes
	// copied is returned.
	uint32 getPixels(uint32 image, uint32 channel, char* dest);
	
	// getVoxelX, Y, and Z, return the corresponding physical dimensions of the side
	// of each voxel in meters.
	float64 getVoxelSizeDx() const {return cz_voxel_size_x;}
	float64 getVoxelSizeDy() const {return cz_voxel_size_y;}
	float64 getVoxelSizeDz() const {return cz_voxel_size_z;}
	
//	class out_of_bounds {}; // For exceptions


//	class Iterator;
	class LSM_Directory;
	typedef std::vector<LSM_Directory> DirectoryContainer;
	typedef DirectoryContainer::iterator iterator;
	// typedef DirectoryContainer::const_iterator const_iterator;

	iterator begin()  {return d_lsmDir.begin();}	
	iterator end() {return d_lsmDir.end();}
	//	const_iterator begin() const {return d_lsmDir.begin();}
	//  const_iterator end() const {return d_lsmDir.end();}

	
	LSM_Directory& operator[](uint32 i) {
		if(i < 0 || i >= d_lsmDir.size()) {
			throw std::out_of_range("LSM_File::operator[i] index is illegal.");
		}
		return d_lsmDir[i];
	}
	
private:	
	
	std::ifstream* d_in;
	uint16 d_lsmByteOrder;
	bool   d_reverseByteOrder;
	uint16 d_magicNumber;
	uint32 d_firstDirectoryOffset;
	//std::vector<LSM_Directory> d_lsmDir;
	DirectoryContainer d_lsmDir;

	// CZ-Private Tags:
	uint32 cz_magic_number;
	sint32 cz_structure_size;
	sint32 cz_dimension_x;
	sint32 cz_dimension_y;
	sint32 cz_dimension_z;
	sint32 cz_dimension_channels;
	sint32 cz_dimension_time;
	sint32 cz_data_type;
	sint32 cz_thumbnail_x;
	sint32 cz_thumbnail_y;
	float64 cz_voxel_size_x;
	float64 cz_voxel_size_y;
	float64 cz_voxel_size_z;
	uint32 cz_scan_type;
	uint32 cz_data_type_2;
	uint32 cz_offset_vector_overlay;
	uint32 cz_offset_input_lut;
	uint32 cz_offset_output_lut;
	uint32 cz_offset_channel_colors;
	float64 cz_time_interval;
	uint32 cz_offset_channel_data_types;
	uint32 cz_offset_scan_information;
	uint32 cz_offset_ks_data;
	uint32 cz_offset_time_stamps;
	uint32 cz_offset_event_list;
	uint32 cz_offset_roi;
	uint32 cz_offset_bleach_roi;
	uint32 cz_offset_next_recording;
	uint32 cz_reserved;
	
	// CZ-Channel Names and Colors Table
	sint32 cnc_block_size;
	sint32 cnc_number_of_colors;
	
	sint32 cnc_number_of_names; // color names, should match cnc_number_of_colors
	sint32 cnc_colors_offset;
	std::vector< uint32 > cnc_colors;
	sint32 cnc_names_offset;
	std::vector< std::string > cnc_color_names;
	sint32 cnc_mono;
	
	// CZ Scan Info
	
};

class LSM_File::LSM_Directory {
public:
	friend class LSM_File;
	
	LSM_Directory();	// Default constructor


	void parseNextTag(char *c, bool reverseByteOrder);
	
	// In TIFF files, some tagged values can represent either values, or a file
	// offset (in the event that the values do not fit in a 32 bit field). The
	// function computeOffsets() parses each of these fields and determines
	// their values.
	void computeOffsets(std::ifstream& in, bool reverseByteOrder);
	
	// Print data members to the indicated output stream.
	void dump(std::ostream &os) const;
	
	// isThumbnail returns true if this image is a low resolution copy
	// (or thumbnail) of another image within the file. A value of false
	// indicates that this image is full resolution.
	bool isThumbnail(void) const;
	
	// isPlanar returns true if the image is in "planar format", in which 
	// each component is stored as a separate plane. A value of false is
	// returned if the components are stored in "chunky format", in which
	// the components are interleaved: RGBRGBRGB ...
	// Note that if there is only one component, there is no distinction
	// between the two. In this case, isPlanar also returns true.
	bool isPlanar(void) const;
	
	// isRGB returns true if this image is an RGB image.
	bool isRGB(void) const;
	
	// getStripByteCount returns the number of bytes in the indicated image
	// channel. 
	uint32 getStripByteCount(uint32 channel) const;
	
	// getStripOffset returns the location of the image data for then indicated
	// channel. 
	uint32 getStripOffset(uint32 channel) const;
	
	// getSamplesPerPixel(void) returns the number of channels in the current
	// image.
	uint16 getSamplesPerPixel(void) const;
	
	// getBitsPerSample returns the number of bits for the indicated channel. 
	// This value is typically 8 for each channel of 24-bit RGB image, for
	// which channel should equal 0, 1, or 2.
	uint16 getBitsPerSample(uint32 channel) const;
	
	// getLength returns the (vertical) length of the current image in pixels.
	uint32 getHeight(void) const;
	uint32 getLength(void) const {return getHeight();}

	// getWidth returns the (horizontal) width of the current image in pixels.
	uint32 getWidth(void) const;
	
	// getPixels copies a block of image data from the current frame in the
	// LSM file, to a block of memory indicated by the second argument. The 
	// number of bytes copied is returned.
    // 
	uint32 getPixels(uint32 channel, char* dest);
	
	// Set the value of the nextDirOffset.
	void setNextDirOffset(uint32 value);
	uint32 getCZOffset(void) const;
	
	// Reset all member data fields to their default (initial) values.
	void reset();
	
private:
	/* Entry keys for the scan information subblock in the LSM_File */
	static const uint32 SUBBLOCK_RECORDING             = 0x10000000;
	static const uint32 SUBBLOCK_LASERS                = 0x30000000;
	static const uint32 SUBBLOCK_LASER                 = 0x50000000;
	static const uint32 SUBBLOCK_TRACKS                = 0x20000000;
	static const uint32 SUBBLOCK_TRACK                 = 0x40000000;
	static const uint32 SUBBLOCK_DETECTION_CHANNELS    = 0x60000000;
	static const uint32 SUBBLOCK_DETECTION_CHANNEL     = 0x70000000;
	static const uint32 SUBBLOCK_ILLUMINATION_CHANNELS = 0x80000000;
	static const uint32 SUBBLOCK_ILLUMINATION_CHANNEL  = 0x90000000;
	static const uint32 SUBBLOCK_BEAM_SPLITTERS        = 0xA0000000;
	static const uint32 SUBBLOCK_BEAM_SPLITTER         = 0xB0000000;
	static const uint32 SUBBLOCK_DATA_CHANNELS         = 0xC0000000;
	static const uint32 SUBBLOCK_DATA_CHANNEL          = 0xD0000000;
	static const uint32 SUBBLOCK_TIMERS                = 0X11000000;
	static const uint32 SUBBLOCK_TIMER                 = 0x12000000;
	static const uint32 SUBBLOCK_MARKERS               = 0x13000000;
	static const uint32 SUBBLOCK_MARKER                = 0x14000000;
	static const uint32 SUBBLOCK_END                   = 0xFFFFFFFF;
	
	static const uint16 LMSTAG_CZ_LSMINFO = 34412;
	
	// This pointer may be redundant since the information of the input stream is maintained
	// by the LSM_File object.
	std::ifstream* d_in;
	
	uint32 d_newSubFileType;
	uint32 d_imageWidth;
	uint32 d_imageHeight;
	
	uint32 d_bitsPerSample_length;  // Also, the number of channels.
	uint32 d_bitsPerSample_offset;
	std::vector<uint16> d_bitsPerSample;
	
	uint16 d_compression;
	uint16 d_photometricinterpretation;
	
	uint32 d_stripoffsets_length;
	uint32 d_stripoffsets_offset;
	std::vector<uint32> d_stripoffsets;
	
	uint16 d_samplesperpixel;
	
	uint32 d_stripbytecounts_length;
	uint32 d_stripbytecounts_offset;
	std::vector<uint32> d_stripbytecounts;
	
	uint32 d_colormap_length;
	uint32 d_colormap_offset;	
	std::vector<uint16> d_colormap;
	
	uint16 d_planarconfiguration;
	uint16 d_predictor;
	uint32 d_cz_lsminfo;
	
	uint32 d_nextDirOffset;
};

// The template function readValue is used to read values of different LSM
// numerical types (i.e., uint8, uint16, uint32, and float64) from either
// a C character string, or an input file stream. Since LSM file data
// is stored in little endian format (the least significant byte is read
// before the most significant byte in multibyte types), and most computers
// assume big endian, these macros facilitate the appropriate byte reversal.

template<typename T>
  T readValue(char *c, const bool reverseByteOrder = false) {
  	if (reverseByteOrder) {
		int i = 0;
		int j = sizeof(T) - 1;
	    for(; i < j; i++, j--) {
	       std::swap(c[i], c[j]);
	    }
	}
	
#if (DEBUG > 99)
	std::cerr << "lsm.h: readValue(char*, ...): value = ";
	for(int i = 0; i < sizeof(T); i++) {
		std::cerr << std::hex << (c[i] & 0xff) << " ";
	}
	std::cerr << std::endl;
#endif

	T* value = reinterpret_cast<T *>(c);  
	return *value;
}

template<typename T>
  T readValue(std::ifstream &in, const bool reverseByteOrder = false) {
	char c[sizeof(T)];
	
	in.read(c, sizeof(T));
	
#if (DEBUG > 99)
	std::cerr << "lsm.h: readValue(ifstream&, ...): reading " 
	          << std::dec << sizeof(T) << " characters: ";
	for(int i = 0; i < sizeof(T); i++) {
		std::cerr << std::hex << (c[i] & 0xff) << " ";
	}
	std::cerr << std::endl;
#endif

	return readValue<T>(c, reverseByteOrder);
}



// Function isLSM_File returns true if and only if the indicated file is an accessible LSM file.
bool isLSM_File(const char* path);

#endif
