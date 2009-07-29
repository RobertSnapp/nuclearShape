/* lsm.cc
 * 
 * Implementation of the classes defined in lsm.h
 *
 * Created on 06 Jan 2007
 *
 */
 
#include <sstream>
#include <stdexcept>
#include "lsm.h"
#include <algorithm> //required for std::swap

// Symbolic constants for photometric interpretation options:
const int TIF_WHITEISZERO = 0;  // Not used by LSM format(?).
const int TIF_BLACKISZERO = 1;
const int TIF_RGB24       = 2;
const int TIF_PALETTE     = 3;

// Symbolic constants for planar configuration options:
const int TIF_CHUNKY = 1;  // default for TIFF 6.0
const int TIF_PLANAR = 2;

using namespace std;

bool lms_type_check() {
	return (sizeof(uint8)   == 1 && 
	        sizeof(uint16)  == 2 &&
	        sizeof(uint32)  == 4 &&
			sizeof(float64) == 8   );
}

// In big-endian format the most significant bytes appear first
bool is_big_endian() {
  unsigned char be_one[] = {0x00, 0x01};
  unsigned short *x = reinterpret_cast<unsigned short *>(be_one);
  return (*x == 1);
}

// In little-endian format the least significant bytes appear first.
bool is_little_endian() {
  unsigned char le_one[] = {0x01, 0x00};
  unsigned short *x = reinterpret_cast<unsigned short *>(le_one);
  return (*x == 1);
}

uint16 getSystemByteOrder() {
// bool isByteOrderReversed() {
// bool reverseByteOrder;
  uint16 byteOrder;

  if (is_little_endian()) {
	// reverseByteOrder = false;
	byteOrder = TIFF_LITTLEENDIAN;
	lsm_debug(std::cerr << "Architecture assumes little_endian byte order.\n");
  } else if (is_big_endian()) {
	// reverseByteOrder = true;
	byteOrder = TIFF_BIGENDIAN;
	lsm_debug(std::cerr << "Architecture assumes big_endian byte order.\n");
  } else {
	lsm_debug(std::cerr << "ERROR: Architecture assume a mixed endian byte order.\n");
	exit(EXIT_FAILURE);
  }
  return byteOrder; // reverseByteOrder;
}


void echoBytes(ifstream &in, int n) {
  unsigned char c;

  std::cerr << "echoBytes:\n";

  
  in.clear();
  for(int i = 0; i < n; i++) {
	in.read((char *) &c, 1);
	if (in.gcount() != 1) {
	  if (in.good()) {
		std::cerr << "Input stream is good.\n";
	  }
	  if (in.bad()) {
		std::cerr << "Input stream is bad.\n";
	  }
	  if (in.eof()) {
		std::cerr << "Input stream is at EOF.";
	  }
	  if (in.eof()) {
		std::cerr << "Input stream failed.";
	  }
	  std::cerr << "Istream failed gcount test: End of file?\n";
	} else {
	// c[0] = in.get();
	  std::cerr << std::hex << (c & 0xff) << " ";
	}
  }
  std::cerr << std::endl;
}

// Default constructor that initializes some members to TIFF 6.0 defaults.
LSM_File::LSM_Directory::LSM_Directory() : d_in(0),
  d_newSubFileType(0), d_imageWidth(0), 
	d_imageHeight(0), d_bitsPerSample_length(0), d_bitsPerSample_offset(0),
	d_compression(0), d_photometricinterpretation(0), d_stripoffsets_length(0),
	d_stripoffsets_offset(0), d_samplesperpixel(1), d_stripbytecounts_length(0),
	d_stripbytecounts_offset(0), d_colormap_length(0), d_colormap_offset(0),
	d_planarconfiguration(1), d_predictor(0), d_cz_lsminfo(0), d_nextDirOffset(0) { }

// Resets the directory fields to their default values:	
void LSM_File::LSM_Directory::reset() {

  lsm_debug(cout << "Calling LSM_File::LSM_Directory::reset()" << endl);

  d_in = 0;
  d_newSubFileType = 0;
  d_imageWidth = 0;
  d_imageHeight = 0;
	
  d_bitsPerSample_length = 0;  // Also, the number of channels.
  d_bitsPerSample_offset = 0;
  d_bitsPerSample.clear();
  d_compression = 0;
  d_photometricinterpretation = 0;
	
  d_stripoffsets_length = 0;
  d_stripoffsets_offset = 0;
  d_stripoffsets.clear();
	
  d_samplesperpixel = 1;
	
  d_stripbytecounts_length = 0;
  d_stripbytecounts_offset = 0;
  d_stripbytecounts.clear();
	
  d_colormap_length = 0;
  d_colormap_offset = 0;	
  d_colormap.clear();
	
  d_planarconfiguration = 0;
  d_predictor = 0;
  d_cz_lsminfo = 0;
	
  d_nextDirOffset = 0;
}

void LSM_File::LSM_Directory::computeOffsets(ifstream& in, bool const byteOrder) {
	if (d_bitsPerSample.size() == 0) {
		if (d_bitsPerSample_length <= 2) {
			d_bitsPerSample.push_back(uint16(d_bitsPerSample_offset >> 16));
			d_bitsPerSample.push_back(uint16(d_bitsPerSample_offset & 0x00ff));
		} else {
			in.seekg(d_bitsPerSample_offset);
			for(uint32 i = 0; i < d_bitsPerSample_length; i++) {
				uint16 bitsPerSample = readValue<uint16>(in, byteOrder);
				d_bitsPerSample.push_back(bitsPerSample);
			}
		}
	} else {
		cerr << "Warning: attempted to overwrite d_bitsPerSample." << endl;
	}
	
	if (d_colormap.size() == 0) {
		if (d_colormap_length <= 2) {
			d_colormap.push_back(uint16(d_colormap_offset >> 16));
			d_colormap.push_back(uint16(d_colormap_offset & 0x00ff));
		} else {
			in.seekg(d_colormap_offset);
			for(uint32 i = 0; i < d_colormap_length; i++) {
				uint16 colorMap = readValue<uint16>(in, byteOrder);
				d_colormap.push_back(colorMap);
			}
		}
	} else {
		cerr << "Warning: attempted to overwrite d_bitsPerSample." << endl;
	}
	
	if (d_stripoffsets.size() == 0) {
		if (d_stripoffsets_length <= 1) {
			d_stripoffsets.push_back(d_stripoffsets_offset);
		} else {
			in.seekg(d_stripoffsets_offset);
			for(uint32 i = 0; i < d_stripoffsets_length; i++) {
				uint32 stripOffsets = readValue<uint32>(in, byteOrder);
				d_stripoffsets.push_back(stripOffsets);
			}
		}
	} else {
		cerr << "Warning: attempted to overwrite d_stripoffsets." << endl;
	}
	
	if (d_stripbytecounts.size() == 0) {
		if (d_stripbytecounts_length <= 1) {
			d_stripbytecounts.push_back(d_stripbytecounts_offset);
		} else {
			in.seekg(d_stripbytecounts_offset);
			for(uint32 i = 0; i < d_stripbytecounts_length; i++) {
				uint32 stripByteCounts = readValue<uint32>(in, byteOrder);
				d_stripbytecounts.push_back(stripByteCounts);
			}
		}
	} else {
		cerr << "Warning: attempted to overwrite d_stripbytecounts." << endl;
	}	
}

void LSM_File::LSM_Directory::dump(ostream &os) const {
	os << "New Subfile = " << dec << d_newSubFileType << endl;
	os << "Image Width = " << dec << d_imageWidth << endl;
	os << "Image Length = " << dec << d_imageHeight << endl;
	
	os << "BitsPerSampleLength = " << d_bitsPerSample_length << endl;
	os << "BitsPerSampleOffset = " << d_bitsPerSample_offset << endl;	
	for(uint32 i = 0; i < d_bitsPerSample.size(); i++) {
		os << "BitsPerSample[" << i << "] = " << d_bitsPerSample[i] << endl;
	}
	
	os << "Compression = " << d_compression << endl;
	os << "Photometric = " << d_photometricinterpretation << endl;
		
	os << "StripOffsetsLength = " << d_stripoffsets_length << endl;
	os << "StripOffsetsOffset = " << hex<< d_stripoffsets_offset << endl;	
	for(uint32 i = 0; i < d_stripoffsets.size(); i++) {
		os << "StripOffsets[" << i << "] = " << hex << d_stripoffsets[i] << endl;
	}
	
	os << "SamplesPerPixel = " << dec << d_samplesperpixel << endl;
	
	os << "StripByteCountsLength = " << dec<< d_stripbytecounts_length << endl;
	os << "StripByteCountssOffset = " << hex << d_stripbytecounts_offset << endl;	
	for(uint32 i = 0; i < d_stripbytecounts.size(); i++) {
		os << "StripByteCounts[" << i << "] = " 
		   << dec << d_stripbytecounts[i] << endl;
	}
	
	os << "ColorMapLength = " << dec << d_colormap_length << endl;
	os << "ColorMapOffset = " << hex << d_colormap_offset << endl;	
	for(uint32 i = 0; i < d_colormap.size(); i++) {
		os << "ColorMap[" << i << "] = " << hex << d_colormap[i] << endl;
	}
	
	os << "PlanarConfiguration = " << hex << d_planarconfiguration << endl;
	os << "Predictor           = " << hex << d_predictor << endl;
	os << "LMS_Info            = " << hex << d_cz_lsminfo << endl;
	os << "NextDirectoryOffset = " << hex << d_nextDirOffset << endl << endl;
}

uint32 LSM_File::LSM_Directory::getCZOffset() const {
	return d_cz_lsminfo;
}

// getBitsPerSample returns the number of bits for the indicated channel. 
// This value is typically 8 for each channel of 24-bit RGB image, for
// which channel should equal 0, 1, or 2.
uint16 LSM_File::LSM_Directory::getBitsPerSample(uint32 channel) const {
	if (channel < 0 || channel >= d_bitsPerSample.size()) {
		throw out_of_range("getBitsPerSample: Illegal channel.");
	}
	return d_bitsPerSample[channel];
}

uint32 LSM_File::LSM_Directory::getStripByteCount(uint32 channel) const {
	if (channel < 0 || channel >= d_stripbytecounts.size()) {
		throw out_of_range("getStripByteCount: Illegal channel.");
	}
	return d_stripbytecounts[channel];	
}

uint32 LSM_File::LSM_Directory::getStripOffset(uint32 channel) const {
	if (channel < 0 || channel >= d_stripoffsets.size()) {
		throw out_of_range("getStripOffset: Illegal channel.");
	}
	return d_stripoffsets[channel];
}

// getSamplesPerPixel(void) returns the number of channels in the current image.
uint16 LSM_File::LSM_Directory::getSamplesPerPixel(void) const {
	return d_samplesperpixel;
}

// isThumbnail returns true if this image is a low resolution copy
// (or thumbnail) of another image within the file. A value of false
// indicates that this image is full resolution.	
bool LSM_File::LSM_Directory::isThumbnail() const {
	return d_newSubFileType & 0x0001;
}

bool LSM_File::LSM_Directory::isPlanar(void) const {
	return (d_planarconfiguration == TIF_PLANAR | d_samplesperpixel == 1);
}
// isRGB returns true if this image is an RGB image.
bool LSM_File::LSM_Directory::isRGB() const {
	return (d_photometricinterpretation == TIF_RGB24);
}

uint32 LSM_File::LSM_Directory::getWidth() const {
	return d_imageWidth;
}

uint32 LSM_File::LSM_Directory::getHeight() const {
	return d_imageHeight;
}

// getPixels copies a block of image data from the current LSM file, to a
// block of memory indicated by the third argument. The number of bytes
// copied is returned.
uint32 LSM_File::LSM_Directory::getPixels(uint32 channel, char* dest) {
	if (channel < 0 || channel >= d_samplesperpixel) 
		throw out_of_range("getPixels: first argument (channel) is illegal.");
	d_in->seekg(getStripOffset(channel));
	uint32 bytesToRead = getStripByteCount(channel);
	(void) d_in->read(dest, bytesToRead);
	return bytesToRead;
}


// Function parseNextTag is called by the LSM_File constructor to parse each
// directory entry.

void LSM_File::LSM_Directory::parseNextTag(char *c, bool const byteOrder) {
	char *ptr = c;
	
	uint16 tag  = readValue<uint16>(ptr, byteOrder);
	ptr += 2;
	uint16 type = readValue<uint16>(ptr, byteOrder);
#pragma unused(type)

	ptr += 2;
	uint32 len  = readValue<uint32>(ptr, byteOrder);
	ptr += 4;
	
	lsm_debug(cerr << "parseNextTag: tag = " << dec << tag
	           << ", type = " << type
	           << ", length = " << len << endl);
	          
	switch (tag) {
		case TIFFTAG_SUBFILETYPE:
			if (len != 1) {
				cerr << "Error: " << len << " is an illegal tag length!" << endl;
				exit(1);
			}
			d_newSubFileType = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_newSubFileType = " << d_newSubFileType << endl);
			break;
			
		case TIFFTAG_IMAGEWIDTH:
			if (len != 1) {
				cerr << "Error: " << len << " is an illegal tag length!" << endl;
				exit(1);
			}
			d_imageWidth = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_imageWidth (dec) = " << dec << d_imageWidth << endl);
			break;
			
		case TIFFTAG_IMAGELENGTH:
			if (len != 1) {
				cerr << "Error: " << len << " is an illegal tag length!" << endl;
				exit(1);
			}
			d_imageHeight = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_imageHeight (dec) = " << dec << d_imageHeight << endl);
			break;
			
		case TIFFTAG_BITSPERSAMPLE:
			d_bitsPerSample_length = len;
			d_bitsPerSample_offset = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_bitsPerSample_length (dec) = " << dec << len << endl);
			break;
			
		case TIFFTAG_COMPRESSION:
			d_compression = readValue<uint16>(ptr, byteOrder);
			if (d_compression != 1) {
				cerr << "Error: This implementation does not "
				     << "support image compression." << endl;
				exit(1);
			}
			lsm_debug(cerr << "d_compression (dec) = " << dec << d_compression << endl);
			break;
			
		case TIFFTAG_PHOTOMETRIC:
			d_photometricinterpretation = readValue<uint16>(ptr, byteOrder);
			lsm_debug(cerr << "d_photometricinterpretation = " 
					   << d_photometricinterpretation << endl);
			break;
			
		case TIFFTAG_STRIPOFFSETS:
			d_stripoffsets_length = len;
			d_stripoffsets_offset = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_stripoffsets_length = " << len << endl);
			break;
			
		case TIFFTAG_SAMPLESPERPIXEL:
			if (len != 1) {
				cerr << "Error: " << len << " is an illegal tag length!" << endl;
				exit(1);
			}
			d_samplesperpixel = readValue<uint16>(ptr, byteOrder);
			lsm_debug(cerr << "d_samplesperpixel = " << d_samplesperpixel << endl);
			break;
			
		case TIFFTAG_STRIPBYTECOUNTS:
			d_stripbytecounts_length = len;
			d_stripbytecounts_offset = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_stripbytecounts_length = " << len << endl);
			break;
			
		case TIFFTAG_COLORMAP:
			d_colormap_length = len;
			d_colormap_offset = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_colormap_length = " << len << endl);
			break;
			
		case TIFFTAG_PLANARCONFIG:
			d_planarconfiguration = readValue<uint16>(ptr, byteOrder);
			lsm_debug(cerr << "d_planarconfiguration = " << d_planarconfiguration << endl);
			break;
			
		case TIFFTAG_PREDICTOR:
			d_predictor = readValue<uint16>(ptr, byteOrder);
			lsm_debug(cerr << "d_predictor (hex) = " << hex << d_predictor << endl);
			break;
			
		case LMSTAG_CZ_LSMINFO:
			d_cz_lsminfo = readValue<uint32>(ptr, byteOrder);
			lsm_debug(cerr << "d_cz_lsminfo (hex) = " << hex << d_cz_lsminfo << endl);
			break;
			
		default:
			cerr << "Error: Unsupported Image Tag = " << tag << endl;
			exit(1);
	}
}

void LSM_File::LSM_Directory::setNextDirOffset(uint32 value) {
	d_nextDirOffset = value;
}


// Function isLSM_File returns true if and only if the indicated file is an accessible LSM file.
bool isLSM_File(const char* path) {
  std::ifstream* input = new ifstream(path, ios::in | ios::binary);
  if (! input) {
    return false;
  }

  uint16 theByteOrder = readValue<uint16>(*input);
  bool reverseByteOrder = (theByteOrder != getSystemByteOrder());

  uint16 theMagicNumber = readValue<uint16>(*input, reverseByteOrder);	
  bool status = (theMagicNumber == LSM_File::LSM_MAGIC_NUMBER);
    
  input->close();
  return status;
}


// constructor of the principal class: LSM_File.
LSM_File::LSM_File(const char* path) : d_in(0), d_lsmByteOrder(0), d_magicNumber(0),
									   d_firstDirectoryOffset(0), cz_magic_number(0),
									   cz_structure_size(0), cz_dimension_x(0),
									   cz_dimension_y(0), cz_dimension_z(0),
									   cz_dimension_channels(0), cz_dimension_time(0),
									   cz_data_type(0), cz_thumbnail_x(0),
									   cz_thumbnail_y(0), cz_voxel_size_x(0),
									   cz_voxel_size_y(0), cz_voxel_size_z(0),
									   cz_scan_type(0), cz_data_type_2(0),
									   cz_offset_vector_overlay(0), cz_offset_input_lut(0),
									   cz_offset_output_lut(0), cz_offset_channel_colors(0),
									   cz_time_interval(0), cz_offset_channel_data_types(0),
									   cz_offset_scan_information(0), cz_offset_ks_data(0),
									   cz_offset_time_stamps(0), cz_offset_event_list(0),
									   cz_offset_roi(0), cz_offset_bleach_roi(0),
									   cz_offset_next_recording(0), cz_reserved(0),
									   cnc_block_size(0), cnc_number_of_colors(0),
									   cnc_number_of_names(0), cnc_colors_offset(0),
									   cnc_names_offset(0), cnc_mono(0)
{  // Opens the LSM file and reads the header.
	if(! lms_type_check()) {
		cerr << "Error: LSM types (uint8, uint16, uint32, float64) "
		     << "have incorrect sizes." << endl;
		exit(EXIT_FAILURE);
	}
	
	d_in = new ifstream(path, ios::in | ios::binary);
	if (! d_in) {
		cerr << "Error: Can't open " << path << endl;
        throw std::ios_base::failure("LSM_File::Can't open fileName.");
	}

	lsm_debug(cerr << "Reading file " << path << endl);

	// Read the lsm byte order from the header. The result should be
	// TIFF_LITTLEENDIAN, a two-byte constant, value that is invariant to order.
	d_lsmByteOrder = readValue<uint16>(*d_in); 
	d_reverseByteOrder = (d_lsmByteOrder != getSystemByteOrder());

	switch (d_lsmByteOrder) {
	case TIFF_LITTLEENDIAN:
	  lsm_debug(cerr << "LSM byte order is little endian." << endl);
	  break;
	case TIFF_BIGENDIAN:
	  lsm_debug(cerr << "LSM byte order is big endian." << endl);
	  break;
	default:
	  cerr << "Illegal value for LSM byte order = " << d_lsmByteOrder << endl;
	  exit(EXIT_FAILURE);
	}

	d_magicNumber = readValue<uint16>(*d_in, d_reverseByteOrder);	
	if (d_magicNumber != LSM_MAGIC_NUMBER) {
		cerr << "Error: File " << path << " is not an LSM file." << endl;
		this->close();
		delete d_in;
        throw std::ios_base::failure("LSM_File: not an LSM file.");
	}	
	lsm_debug(cerr << "File is an LSM file." << endl);
		
	uint32 nextDirOffset = d_firstDirectoryOffset 
		= readValue<uint32>(*d_in, d_reverseByteOrder);
	
	int dircount = 0;
	while (nextDirOffset) {
		lsm_debug(cerr << dec << dircount << " directory offset = " 
		               << hex << nextDirOffset << endl);
		LSM_Directory nextdir; // The working lsm directory structure that is used to build
		                       // incrementally, the d_lsmDir vector.

		// Reset nextdir data to default values.
		if (dircount > 0) nextdir.reset();

		nextdir.d_in = d_in;
		
		d_in->seekg(nextDirOffset);

		// Each lsm image directory begins with a uint16 number that specifies the number
		// of 12 byte entries that immediately follow. Each 12 byte directory entry assumes
		// the following form:
		//     1. A two byte tag that defines the quantity (uint16),
		//     2. A two byte type index (uint16): 
		//        TIF_BYTE = 0x0001, TIF_ASCII = 0x0002, TIF_SHORT = 0x0003, TIF_LONG = 0x0004,
		//        and TIF_RATIONAL = 0x0005.
		//     3. A four byte length (uint32),
		//     4. A four byte value (uint32). In the event that the value requires more than
		//        four bytes, then this field contains a four-byte offset to a location, from
        //        which the complete value can be read.

		uint16 ntags = readValue<uint16>(*d_in, d_reverseByteOrder);
		lsm_debug(cerr << "Reading " << dec << ntags << "tags:" << endl);
		for (int i = 0; i < ntags; i++) {
			char c[12];
			d_in->read(c, 12);
			
			for(int j = 0; j < 12; j++) {
				lsm_debug(cerr << hex << static_cast<int>(c[j] & 0xff) << " ");
			}
			lsm_debug(cerr << endl);

			// parse the 12 byte directory entry.
			nextdir.parseNextTag(c, d_reverseByteOrder);		
		}

		// The last directory item contains a uint32 offset to the next image directory.
		// The last item in the last directory should indicate an offset of 0.
		nextDirOffset = readValue<uint32>(*d_in, d_reverseByteOrder);
		nextdir.setNextDirOffset(nextDirOffset);
		
		nextdir.computeOffsets(*d_in, d_reverseByteOrder);
		if (dircount == 0) {
			uint32 cz_offset = nextdir.getCZOffset();
			parseCZInfo(cz_offset);
			if (cz_offset_channel_colors) {
				parseColorTable(cz_offset_channel_colors);
			}	
		}		
		d_lsmDir.push_back(nextdir);
		lsm_debug(nextdir.dump(cerr));
		dircount++;	
	}	
}

// Destructor: closes the file stream.
LSM_File::~LSM_File() {	
  lsm_debug(cout << "Calling LSM_File::~LSM_File()" << endl);

  this->close();
  if (d_in) {
	delete d_in;
	d_in = 0;
  }
  d_lsmDir.clear();
}

// getImageCount returns the number of images found in the current file.
uint32 LSM_File::getImageCount(void) const {
	return d_lsmDir.size();
}

// isThumbnail returns true if the indicated image is a low 
// resolution copy (or thumbnail) of another image within the file. 
// A value of false indicates that this image is full resolution.
bool LSM_File::isThumbnail(uint32 image) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::isThumbnail: Illegal image index.");
	}
	return d_lsmDir[image].isThumbnail();
}

// isPlanar returns true if the indicated image is in "planar format",
// in which each component is stored as a separate plane. A value of 
// false is returned if the components are stored in "chunky format",
// in which the components are interleaved: RGBRGBRGB ...
// Note that if there is only one component, there is no distinction
// between the two. In this case, isPlanar also returns true.
bool LSM_File::isPlanar(uint32 image) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::isPlanar: Illegal image index.");
	}
	return d_lsmDir[image].isPlanar();
}

// isRGB returns true if the indicated image is stored in RGB format.
bool LSM_File::isRGB(uint32 image) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::isRGB: Illegal image index.");
	}
	return d_lsmDir[image].isRGB();
}

// getStripByteCount returns the number of bytes in the indicated image
// channel for the indicated image.
uint32 LSM_File::getStripByteCount(uint32 image, uint32 channel) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::getStripByteCount: Illegal image index.");
	}
	if (channel < 0 || channel >= getSamplesPerPixel(image)) {
		throw out_of_range("LSM_File::getStripByteCount: Illegal channel.");
	}
	return d_lsmDir[image].getStripByteCount(channel);
}

// getStripOffset returns the location of the image data for then indicated
// channel in the indicated image.
uint32 LSM_File::getStripOffset(uint32 image, uint32 channel) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::getStripOffset: Illegal image index.");
	}
	if (channel < 0 || channel >= getSamplesPerPixel(image)) {
		throw out_of_range("LSM_File::getStripOffset: Illegal channel.");
	}
	return d_lsmDir[image].getStripOffset(channel);
}

// getSamplesPerPixel(void) returns the number of channels in the
// indicated image.
uint16 LSM_File::getSamplesPerPixel(uint32 image) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::getSamplesPerPixel: Illegal image index.");
	}
	return d_lsmDir[image].getSamplesPerPixel();
}

// getBitsPerSample returns the number of bits for the indicated channel
// in the indicated image. This value is typically 8 for each channel of
// 24-bit RGB image, for which channel should equal 0, 1, or 2.
uint16 LSM_File::getBitsPerSample(uint32 image, uint32 channel) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::getBitsPerSample: Illegal image index.");
	}
	if (channel < 0 || channel >= getSamplesPerPixel(image)) {
		throw out_of_range("LSM_File::getBitsPerSample: Illegal channel.");
	}
	return d_lsmDir[image].getBitsPerSample(channel);
}



// getHeight returns the (vertical) length of the indicated image in pixels.
uint32 LSM_File::getHeight(uint32 image) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::getHeight: Illegal image index.");
	}
	return d_lsmDir[image].getHeight();
}

// getWidth returns the (horizontal) width of the indicated image in pixels.
uint32 LSM_File::getWidth(uint32 image) const {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("LSM_File::getWidth: Illegal image index.");
	}
	return d_lsmDir[image].getWidth();
}

// getPixels copies a block of image data from the current LMS file, to a
// block of memory indicated by the third argument. The number of bytes
// copied is returned.
uint32 LSM_File::getPixels(uint32 image, uint32 channel, char* dest) {
	if (image < 0 || image >= d_lsmDir.size()) {
		throw out_of_range("getPixels: Illegal image index.");
	}
	if (channel < 0 || channel >= getSamplesPerPixel(image)) {
	  out_of_range("getPixels: Illegal channel.");
	}
	d_in->seekg(d_lsmDir[image].getStripOffset(channel));
	uint32 bytesToRead = d_lsmDir[image].getStripByteCount(channel);
	(void) d_in->read(dest, bytesToRead);
	return bytesToRead;
}


void LSM_File::parseColorTable(uint32 offset) {
	d_in->seekg(offset);
	cnc_block_size       = readValue<uint32>(*d_in, d_reverseByteOrder);
	cnc_number_of_colors = readValue<uint32>(*d_in, d_reverseByteOrder);
	cnc_number_of_names  = readValue<uint32>(*d_in, d_reverseByteOrder);
	cnc_colors_offset    = readValue<uint32>(*d_in, d_reverseByteOrder);
	cnc_names_offset     = readValue<uint32>(*d_in, d_reverseByteOrder);
	cnc_mono             = readValue<uint32>(*d_in, d_reverseByteOrder);
	
	if (cnc_colors_offset) {
		d_in->seekg(cnc_colors_offset);
		for(int i = 0; i < cnc_number_of_colors; i++) {
			uint32 color = readValue<uint32>(*d_in, d_reverseByteOrder);
			cnc_colors.push_back(color);
		}
	}
	
	if (cnc_names_offset) {
		d_in->seekg(cnc_names_offset);
		for(int i = 0; i < cnc_number_of_names; i++) {
			string name;
			char c;
			do {
				d_in->get(c);
				name.push_back(c);
			} while (c);
			cnc_color_names.push_back(name);
			name.clear();
		}
	}
}

void LSM_File::parseCZInfo(uint32 offset) {
	if (! offset) {
		cerr << "Warning: No CZ information is available";
		return;
	}
	
	d_in->seekg(offset);
	cz_magic_number = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "Magic number = " << hex << cz_magic_number << endl);
	
	cz_structure_size = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "Structure Size = " << dec << cz_structure_size << endl);
	
	cz_dimension_x = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "dimension_x = " << dec << cz_dimension_x << endl);
	
	cz_dimension_y = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "dimension_y = " << dec << cz_dimension_y << endl);
	
	cz_dimension_z = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "dimension_z = " << dec << cz_dimension_z << endl);
	
	cz_dimension_channels = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "dimension_channels = " << dec << cz_dimension_channels << endl);
	
	cz_dimension_time = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "dimension_time= " << dec << cz_dimension_time << endl);
	
	cz_data_type = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "data_type = " << dec << cz_data_type << endl);
	
	cz_thumbnail_x = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "thumbnail_x = " << dec << cz_thumbnail_x << endl);
		
	cz_thumbnail_y = readValue<uint32>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "thumbnail_y = " << dec << cz_thumbnail_y << endl);
	
	cz_voxel_size_x = readValue<float64>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "voxel_size_x = " << dec << cz_voxel_size_x << endl);
	
	cz_voxel_size_y = readValue<float64>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "voxel_size_y = " << dec << cz_voxel_size_y << endl);
	
	cz_voxel_size_z = readValue<float64>(*d_in, d_reverseByteOrder);
	lsm_debug(cerr << "voxel_size_z = " << dec << cz_voxel_size_z << endl);
	
	cz_scan_type                 = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_data_type_2               = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_vector_overlay     = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_input_lut          = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_output_lut         = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_channel_colors     = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_time_interval             = readValue<float64>(*d_in, d_reverseByteOrder);
	cz_offset_channel_data_types = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_scan_information   = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_ks_data            = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_time_stamps        = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_event_list         = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_roi                = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_bleach_roi         = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_offset_next_recording     = readValue<uint32>(*d_in, d_reverseByteOrder);
	cz_reserved                  = readValue<uint32>(*d_in, d_reverseByteOrder);
}
