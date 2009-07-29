/*****************************************************************************/
/* pngImage.h 
/* 
/* Defines a class for reading and manipulating a png image.
/*
/* Robert R. Snapp  Copyright (C) 2006.
/*****************************************************************************/

#ifndef __PNG_GRAYIMAGE_H__
#define __PNG_GRAYIMAGE_H__

#include <algorithm> // min function
#include <cmath> // log function
#include <valarray>
#include <GLUT/glut.h>
#include "png.h"
#include <cstdio> // libpng routnes require a c-style file pointer.
#include <iostream>
#include <string>
#include <stdlib.h>
#include <csetjmp>
#include <cassert>
#include <typeinfo>
// #include "pixelArray.h"
#include <fstream>
#include "grayImage.h"

typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;

template<typename T> 
class pngGrayImage : public grayImage<T> {
 protected:
  png_structp png_ptr_;
  png_infop   info_ptr_;
  int pngBitDepth_;
  int pngColorType_;
 public:
 pngGrayImage(int rows, int cols) : grayImage<T>(rows, cols) {
	png_ptr_  = 0;
	info_ptr_ = 0;
  }

 pngGrayImage() : grayImage<T>() {
	png_ptr_ = 0;
	info_ptr_ = 0;
  }
	
  ~pngGrayImage() {}
		
 pngGrayImage(const pngGrayImage &p) : grayImage<T>(grayImage<T>::p) {		
	png_ptr_ = p.png_ptr_;
	info_ptr_ = p.info_ptr_;
  }
		
  void clear() {
	grayImage<T>::clear();
	if (png_ptr_) png_destroy_read_struct(&png_ptr_, &info_ptr_, 0);
  }
	
  pngGrayImage operator=(const pngGrayImage &p) {
	if (this != &p) {
	  png_ptr_ = p.png_ptr_;
	  info_ptr_ = p.info_ptr_;
	}
	return *this;
  }
	
  bool isPNG_File(const char* filename);
  bool readPNGfile(const char* filename);
  void writePNGfile(const char* filename);
 protected:
  void adjust_gamma(double display_exponent);
  int readpng_init(FILE *infile);
  int readpng_get_bgcolor(uchar *red, uchar *green, uchar *blue);
};

/*---------------------------------------------------------------------------

   rpng - simple PNG display program                              readpng.c

  ---------------------------------------------------------------------------

      Copyright (c) 1998-2000 Greg Roelofs.  All rights reserved.

      This software is provided "as is," without warranty of any kind,
      express or implied.  In no event shall the author or contributors
      be held liable for any damages arising in any way from the use of
      this software.

      Permission is granted to anyone to use this software for any purpose,
      including commercial applications, and to alter it and redistribute
      it freely, subject to the following restrictions:

      1. Redistributions of source code must retain the above copyright
         notice, disclaimer, and this list of conditions.
      2. Redistributions in binary form must reproduce the above copyright
         notice, disclaimer, and this list of conditions in the documenta-
         tion and/or other materials provided with the distribution.
      3. All advertising materials mentioning features or use of this
         software must display the following acknowledgment:

            This product includes software developed by Greg Roelofs
            and contributors for the book, "PNG: The Definitive Guide,"
            published by O'Reilly and Associates.

  ---------------------------------------------------------------------------*/



/* future versions of libpng will provide this macro: */
#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr_)   ((png_ptr_)->jmpbuf)
#endif

// void readpng_version_info(void)
// {
//     fprintf(stderr, "   Compiled with libpng %s; using libpng %s.\n",
//       PNG_LIBPNG_VER_STRING, png_libpng_ver);
//     fprintf(stderr, "   Compiled with zlib %s; using zlib %s.\n",
//       ZLIB_VERSION, zlib_version);
// }

    /* unlike the example in the libpng documentation, we have *no* idea where
     * this file may have come from--so if it doesn't have a file gamma, don't
     * do any correction ("do no harm") */
template<typename T>
void pngGrayImage<T>::adjust_gamma(double display_exponent) {
	double gamma;
	
    if (png_get_gAMA(png_ptr_, info_ptr_, &gamma))
        png_set_gamma(png_ptr_, display_exponent, gamma);
}

/* return value = 0 for success, 1 for bad sig, 2 for bad IHDR, 4 for no mem */

template<typename T>
int pngGrayImage<T>::readpng_init(FILE *infile)
{
    uchar sig[8];
	uint cols, rows, pngBitDepth, pngColorType

    /* first do a quick check that the file really is a PNG image; could
     * have used slightly more general png_sig_cmp() function instead */
    fread(sig, 1, 8, infile);
    if (!png_check_sig(sig, 8))
        return 1;   /* bad signature */


    /* could pass pointers to user-defined error handlers instead of NULLs: */
    png_ptr_ = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    if (!png_ptr_)
        return 4;   /* out of memory */

    info_ptr_ = png_create_info_struct(png_ptr_);
    if (!info_ptr_) {
        png_destroy_read_struct(&png_ptr_, 0, 0);
        return 4;   /* out of memory */
    }

    /* we could create a second info struct here (end_info), but it's only
     * useful if we want to keep pre- and post-IDAT chunk info separated
     * (mainly for PNG-aware image editors and converters) */

    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */
    if (setjmp(png_jmpbuf(png_ptr_))) {
        png_destroy_read_struct(&png_ptr_, &info_ptr_, 0);
        return 2;
    }

    png_init_io(png_ptr_, infile);
    png_set_sig_bytes(png_ptr_, 8);  /* we already read the 8 signature bytes */
    png_read_info(png_ptr_, info_ptr_);  /* read all PNG info up to image data */

    /* alternatively, could make separate calls to png_get_image_width(),
     * etc., but want pngBitDepth_ and pngColorType_ for later [don't care about
     * compression_type and filter_type => NULLs] */
    png_get_IHDR(png_ptr_, info_ptr_, &cols_, &rows_, &pngBitDepth_, &pngColorType_,
      0, 0, 0);
	
    /* OK, that's all we need for now; return happy */

    return 0;
}

/* returns 0 if succeeds, 1 if fails due to no bKGD chunk, 2 if libpng error;
 * scales values to 8-bit if necessary */

template<typename T>
int pngGrayImage<T>::readpng_get_bgcolor(uchar *red, uchar *green, uchar *blue)
{
    png_color_16p pBackground;


    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */

    if (setjmp(png_jmpbuf(png_ptr_))) {
        png_destroy_read_struct(&png_ptr_, &info_ptr_, NULL);
        return 2;
    }


    if (!png_get_valid(png_ptr_, info_ptr_, PNG_INFO_bKGD))
        return 1;

    /* it is not obvious from the libpng documentation, but this function
     * takes a pointer to a pointer, and it always returns valid red, green
     * and blue values, regardless of color_type: */

    png_get_bKGD(png_ptr_, info_ptr_, &pBackground);


    /* however, it always returns the raw bKGD data, regardless of any
     * bit-depth transformations, so check depth and adjust if necessary */

    if (pngBitDepth_ == 16) {
        *red   = pBackground->red   >> 8;
        *green = pBackground->green >> 8;
        *blue  = pBackground->blue  >> 8;
    } else if (pngColorType_ == PNG_COLOR_TYPE_GRAY && pngBitDepth_ < 8) {
        if (pngBitDepth_ == 1)
            *red = *green = *blue = pBackground->gray? 255 : 0;
        else if (pngBitDepth_ == 2)
            *red = *green = *blue = (255/3) * pBackground->gray;
        else /* pngBitDepth_ == 4 */
            *red = *green = *blue = (255/15) * pBackground->gray;
    } else {
        *red   = (uchar)pBackground->red;
        *green = (uchar)pBackground->green;
        *blue  = (uchar)pBackground->blue;
    }

    return 0;
}


template<typename T>
bool pngGrayImage<T>::readPNGfile(const char* filename) {
	FILE *fp;
	png_uint_32 i;
	png_uint_32 bands;

	
	fp = fopen(filename, "rb");
	if (! fp) {
		fprintf(stderr, "Can't open file %s\n", filename);
		exit(1);
	}
	
	if (readpng_init(fp) != 0) {
	  // error: file is either not a PNG file, or is corrupted
	  return false;
	}
  
    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */

    if (setjmp(png_jmpbuf(png_ptr_))) {
        png_destroy_read_struct(&png_ptr_, &info_ptr_, 0);
        fprintf(stderr, "readPNGfile: setjmp error");
		return false;
    }


    /* DON'T expand palette images to RGB, low-bit-depth grayscale images to 8 bits,
     * transparency chunks to full alpha channel; strip 16-bit-per-sample
     * images to 8 bits per sample; and convert grayscale to RGB[A] */

    if (pngColorType_ == PNG_COLOR_TYPE_PALETTE)
        png_set_expand(png_ptr_);
    if (pngColorType_ == PNG_COLOR_TYPE_GRAY && pngBitDepth_ < 8)
        png_set_expand(png_ptr_);
    if (png_get_valid(png_ptr_, info_ptr_, PNG_INFO_tRNS))
        png_set_expand(png_ptr_);
    if (pngBitDepth_ == 16)
        png_set_strip_16(png_ptr_);
    if (pngColorType_ == PNG_COLOR_TYPE_GRAY ||
        pngColorType_ == PNG_COLOR_TYPE_GRAY_ALPHA)
        //png_set_gray_to_rgb(png_ptr_);

    /* unlike the example in the libpng documentation, we have *no* idea where
     * this file may have come from--so if it doesn't have a file gamma, don't
     * do any correction ("do no harm") */

	adjust_gamma(1.0);

    /* all transformations have been registered; now update info_ptr data,
     * get bytesPerRow_ and channels, and allocate image memory */

    png_read_update_info(png_ptr_, info_ptr_);

    bytesPerRow_ = png_get_rowbytes(png_ptr_, info_ptr_);   
    bands = (int)png_get_channels(png_ptr_, info_ptr_);
    if (bands != 1) {
    	std::cerr << "grayImage::readPNG: (WARNING) band mismatch." << std::endl;
    }
	cols_ = bytesPerRow_/bands;
	pixel_.resize(bytesPerRow_*rows_);
	
	png_bytep* row_pointers = new png_bytep [rows_];  
    for (i = 0; i < rows_; i++)
        row_pointers[i] = &pixel_[0] + i*bytesPerRow_;

    /* now we can go ahead and just read the whole image */

    png_read_image(png_ptr_, row_pointers);
	delete [] row_pointers;

    png_read_end(png_ptr_, 0);
    fclose(fp);
	return true;
}

template<typename T>
void pngGrayImage<T>::writePNGfile(const char* filename) {
	FILE *fp;
	const int bitsPerByte = 8;
	
	fp = fopen(filename, "wb");
	if (! fp) {
		fprintf(stderr, "pngGrayImage<T>::writePNGfile: Cannot write to file %s in binary mode.\n", filename);
		exit(1);
	}
	
	/* initialize the png file structures (PNG_LIBPNG_VER_STRING is defined
	 * in png.h).
	 */
	png_ptr_ = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      						(void*)         0,  // (void *)user_error_ptr, 
      						(png_error_ptr) 0,  // user_error_fn, 
      						(png_error_ptr) 0); // user_warning_fn);
      						
  if (! png_ptr_) {
  	fclose(fp);
  	fprintf(stderr, "writePNGfile: Cannot create a png write structure.\n");
  	exit (1);
  }

  info_ptr_ = png_create_info_struct(png_ptr_);
  if (! info_ptr_) {
  	fclose(fp);
  	png_destroy_write_struct(&png_ptr_, (png_infopp) NULL);
  	fprintf(stderr, "writePNGfile: Cannot create png info structure.\n");
  	exit(1);
  }

  
  /* libpng handles errors via calls to setjmp and longjmp. Here we create
   * the user error handler with a call to setjmp.
   */ 
   if (setjmp(png_ptr_->jmpbuf)) {
   		/* Activated only if libpng calls longjmp in the future */
   		png_destroy_write_struct(&png_ptr_, &info_ptr_);
   		fprintf(stderr, "writePNGfile: Error encountered in libpng.\n");
  		fclose(fp);
  		exit(1);
  }
 
  /* Tell libpng where to write the png file */
  png_init_io(png_ptr_, fp);
  
  /* Turn the filtering off */
  png_set_filter(png_ptr_, 0, PNG_FILTER_NONE);
  
  /* Prepare the header info */
  info_ptr_->width         = cols_;
  info_ptr_->height        = rows_;
  info_ptr_->bit_depth     = bitsPerByte*sizeof(T);
  info_ptr_->sig_bit.gray   = info_ptr_->bit_depth;
  info_ptr_->color_type    = PNG_COLOR_TYPE_GRAY;
  info_ptr_->interlace_type = 0; /* 0 => None, 1 => Interlaced */
  
  /* Write the header */
  png_write_info(png_ptr_, info_ptr_);
  
  png_bytep* row_pointers = new png_bytep [rows_];
  for (int i = 0; i < rows_; i++)
      row_pointers[i] = &pixel_[0] + i*bytesPerRow_;
        
  /* Process the image data */
  png_write_image(png_ptr_, row_pointers);
  delete [] row_pointers;
  
  png_write_end(png_ptr_, info_ptr_);
  
  /* Clean up */
  	png_destroy_write_struct(&png_ptr_, &info_ptr_);
  	png_ptr_ = 0;
  	info_ptr_ = 0;
	fclose(fp);
}

#endif /* __PNGIMAGE_H__ */