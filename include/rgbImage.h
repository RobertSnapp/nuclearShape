/*****************************************************************************/
/* rgbImage.h 
/* 
/* Robert R. Snapp Copyright (C) 2006.
/*****************************************************************************/

#ifndef __RGBIMAGE_H__
#define __RGBIMAGE_H__

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
#include "pixelArray.h"


typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;


const int RGB_bands = 3;

template<typename T> 
class rgbBase : public pixelArray {
protected:
	std::valarray<T> d_pixel;
	png_structp d_pngPtr;
	png_infop   d_infoPtr;
	int d_pngBitDepth;
	int d_pngColorType;
public:
	rgbBase(int rows, int cols) : pixelArray(rows, cols, RGB_bands) {
		assert(rows > 0);
		assert(cols > 0);
		pixelArray::d_bytesPerRow = d_cols*RGB_bands*sizeof(T);		
		d_pixel.resize(d_rows*d_cols*d_bands);		
		d_pngPtr = 0;
		d_infoPtr = 0;
	}

	rgbBase() : pixelArray() {
		pixelArray::d_bands = RGB_bands;
		d_bytesPerRow = 0;
		d_pixel = 0;
		d_pngPtr = 0;
		d_infoPtr = 0;
	}
	
	~rgbBase() {}
		
	rgbBase(rgbBase & p) {		
		d_pixel.resize(d_rows*d_cols*RGB_bands);
		d_pixel = p.d_pixel;  // Deep copy.
		d_bytesPerRow = d_cols*RGB_bands*sizeof(T);
		d_pngPtr = p.d_pngPtr;
		d_infoPtr = p.d_infoPtr;
	}
		
	rgbBase operator=(rgbBase &p) {
		d_pixel.resize(d_rows*d_cols*RGB_bands);
		d_pixel = p.d_pixel;
		d_pngPtr = p.d_pngPtr;
		d_infoPtr = p.d_infoPtr;
	}
	
	void clear() {
		d_pixel.resize(0);
		if (d_pngPtr) png_destroy_read_struct(&d_pngPtr, &d_infoPtr, 0);
		pixelArray::clear();
	}
	
	T& red(int i, int j) {
		return d_pixel[(i*d_cols + j)*RGB_bands];
	}
	
	T& green(int i, int j) {
		return d_pixel[(i*d_cols + j)*RGB_bands + 1];
	}
	
	T& blue(int i, int j) {
		return d_pixel[(i*d_cols + j)*RGB_bands + 2];
	}
		
    void setPixel(int i, int j, std::valarray<T> &v) {
    	d_pixel[std::slice((i*d_cols + j)*RGB_bands, RGB_bands, 1)] = v;
    }
    
    void getPixel(int i, int j, std::valarray<T> &v) {
    	v = d_pixel[std::slice((i*d_cols + j)*RGB_bands, RGB_bands, 1)];
    }
    
     
    void setRows(int rows) {
    	assert(rows > 0); 
    	d_rows = rows;
    }
    void setCols(int cols) {
    	assert(cols > 0); 
    	d_cols = cols;
    }
   
    void resize() {
    	d_bands = RGB_bands;
    	d_pixel.resize(d_rows*d_cols*RGB_bands);
		d_bytesPerRow = d_cols*RGB_bands*sizeof(T);
	}
    	
    virtual void drawPixels() = 0; 
    	
    void readPNGfile(const char* filename);
    void writePNGfile(const char* filename);
protected:
	void adjust_gamma(double display_exponent);
	int readpng_init(FILE *infile);
	int readpng_get_bgcolor(uchar *red, uchar *green, uchar *blue);

};

template<typename T>
class rgbImage : public rgbBase<T> {
	//rgbImage<T>() : rgbBase<T>() {}
	//rgbImage<T>(int r, int c) : rgbBase<T>(r, c) {}
};

template<>
class rgbImage<GLubyte> : public rgbBase<GLubyte> {
public:
	//rgbImage<GLubyte>() : rgbBase<GLubyte>() {}
	//rgbImage<GLubyte>(int r, int c) : rgbBase<GLubyte>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_UNSIGNED_BYTE, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLshort> : public rgbBase<GLshort> {
public:
	//rgbImage<GLshort>() : rgbBase<GLshort>() {}
	//rgbImage<GLshort>(int r, int c) : rgbBase<GLshort>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_SHORT, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLint> : public rgbBase<GLint> {
public:
	//rgbImage<GLint>() : rgbBase<GLint>() {}
	//rgbImage<GLint>(int r, int c) : rgbBase<GLint>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_INT, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLfloat> : public rgbBase<GLfloat> {
public:
	//rgbImage<GLfloat>() : rgbBase<GLfloat>() {}
	//rgbImage<GLfloat>(int r, int c) : rgbBase<GLfloat>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_FLOAT, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLdouble> : public rgbBase<GLdouble> {
public:
	//rgbImage<GLdouble>() : rgbBase<GLdouble>() {}
	//rgbImage<GLdouble>(int r, int c) : rgbBase<GLdouble>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_DOUBLE, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLbyte> : public rgbBase<GLbyte> {
public:
	//rgbImage<GLbyte>() : rgbBase<GLbyte>() {}
	//rgbImage<GLbyte>(int r, int c) : rgbBase<GLbyte>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_BYTE, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLuint> : public rgbBase<GLuint> {
public:
	//rgbImage<GLuint>() : rgbBase<GLuint>() {}
	//rgbImage<GLuint>(int r, int c) : rgbBase<GLuint>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_UNSIGNED_INT, &d_pixel[0]);
	}
};

template<>
class rgbImage<GLushort> : public rgbBase<GLushort> {
public:
	//rgbImage<GLushort>() : rgbBase<GLushort>() {}
	//rgbImage<GLushort>(int r, int c) : rgbBase<GLushort>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_RGB, GL_UNSIGNED_SHORT, &d_pixel[0]);
	}
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
#  define png_jmpbuf(d_pngPtr)   ((d_pngPtr)->jmpbuf)
#endif

void readpng_version_info(void)
{
    fprintf(stderr, "   Compiled with libpng %s; using libpng %s.\n",
      PNG_LIBPNG_VER_STRING, png_libpng_ver);
    fprintf(stderr, "   Compiled with zlib %s; using zlib %s.\n",
      ZLIB_VERSION, zlib_version);
}

    /* unlike the example in the libpng documentation, we have *no* idea where
     * this file may have come from--so if it doesn't have a file gamma, don't
     * do any correction ("do no harm") */
template<typename T>
void rgbBase<T>::adjust_gamma(double display_exponent) {
	double gamma;
	
    if (png_get_gAMA(d_pngPtr, d_infoPtr, &gamma))
        png_set_gamma(d_pngPtr, display_exponent, gamma);
}

/* return value = 0 for success, 1 for bad sig, 2 for bad IHDR, 4 for no mem */

template<typename T>
int rgbBase<T>::readpng_init(FILE *infile)
{
    uchar sig[8];


    /* first do a quick check that the file really is a PNG image; could
     * have used slightly more general png_sig_cmp() function instead */

    fread(sig, 1, 8, infile);
    if (!png_check_sig(sig, 8))
        return 1;   /* bad signature */


    /* could pass pointers to user-defined error handlers instead of NULLs: */
    d_pngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    if (!d_pngPtr)
        return 4;   /* out of memory */

    d_infoPtr = png_create_info_struct(d_pngPtr);
    if (!d_infoPtr) {
        png_destroy_read_struct(&d_pngPtr, 0, 0);
        return 4;   /* out of memory */
    }


    /* we could create a second info struct here (end_info), but it's only
     * useful if we want to keep pre- and post-IDAT chunk info separated
     * (mainly for PNG-aware image editors and converters) */


    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */

    if (setjmp(png_jmpbuf(d_pngPtr))) {
        png_destroy_read_struct(&d_pngPtr, &d_infoPtr, 0);
        return 2;
    }

    png_init_io(d_pngPtr, infile);
    png_set_sig_bytes(d_pngPtr, 8);  /* we already read the 8 signature bytes */

    png_read_info(d_pngPtr, d_infoPtr);  /* read all PNG info up to image data */


    /* alternatively, could make separate calls to png_get_image_width(),
     * etc., but want d_pngBitDepth and d_pngColorType for later [don't care about
     * compression_type and filter_type => NULLs] */

    png_get_IHDR(d_pngPtr, d_infoPtr, &d_cols, &d_rows, &d_pngBitDepth, &d_pngColorType,
      0, 0, 0);
  


    /* OK, that's all we need for now; return happy */

    return 0;
}

/* returns 0 if succeeds, 1 if fails due to no bKGD chunk, 2 if libpng error;
 * scales values to 8-bit if necessary */

template<typename T>
int rgbBase<T>::readpng_get_bgcolor(uchar *red, uchar *green, uchar *blue)
{
    png_color_16p pBackground;


    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */

    if (setjmp(png_jmpbuf(d_pngPtr))) {
        png_destroy_read_struct(&d_pngPtr, &d_infoPtr, 0);
        return 2;
    }


    if (!png_get_valid(d_pngPtr, d_infoPtr, PNG_INFO_bKGD))
        return 1;

    /* it is not obvious from the libpng documentation, but this function
     * takes a pointer to a pointer, and it always returns valid red, green
     * and blue values, regardless of color_type: */

    png_get_bKGD(d_pngPtr, d_infoPtr, &pBackground);


    /* however, it always returns the raw bKGD data, regardless of any
     * bit-depth transformations, so check depth and adjust if necessary */

    if (d_pngBitDepth == 16) {
        *red   = pBackground->red   >> 8;
        *green = pBackground->green >> 8;
        *blue  = pBackground->blue  >> 8;
    } else if (d_pngColorType == PNG_COLOR_TYPE_GRAY && d_pngBitDepth < 8) {
        if (d_pngBitDepth == 1)
            *red = *green = *blue = pBackground->gray? 255 : 0;
        else if (d_pngBitDepth == 2)
            *red = *green = *blue = (255/3) * pBackground->gray;
        else /* d_pngBitDepth == 4 */
            *red = *green = *blue = (255/15) * pBackground->gray;
    } else {
        *red   = (uchar)pBackground->red;
        *green = (uchar)pBackground->green;
        *blue  = (uchar)pBackground->blue;
    }

    return 0;
}


template<typename T>
void rgbBase<T>::readPNGfile(const char* filename) {
	FILE *fp;
	png_uint_32 i;
	png_uint_32 bands;

	
	fp = fopen(filename, "rb");
	if (! fp) {
		fprintf(stderr, "Can't open file %s\n", filename);
		exit(1);
	}
	
	readpng_init(fp);
  
    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */

    if (setjmp(png_jmpbuf(d_pngPtr))) {
        png_destroy_read_struct(&d_pngPtr, &d_infoPtr, 0);
        fprintf(stderr, "readPNGfile: setjmp error");
       	exit(1);
    }


    /* expand palette images to RGB, low-bit-depth grayscale images to 8 bits,
     * transparency chunks to full alpha channel; strip 16-bit-per-sample
     * images to 8 bits per sample; and convert grayscale to RGB[A] */

    if (d_pngColorType == PNG_COLOR_TYPE_PALETTE)
        png_set_expand(d_pngPtr);
    if (d_pngColorType == PNG_COLOR_TYPE_GRAY && d_pngBitDepth < 8)
        png_set_expand(d_pngPtr);
    if (png_get_valid(d_pngPtr, d_infoPtr, PNG_INFO_tRNS))
        png_set_expand(d_pngPtr);
    if (d_pngBitDepth == 16)
        png_set_strip_16(d_pngPtr);
    if (d_pngColorType == PNG_COLOR_TYPE_GRAY ||
        d_pngColorType == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(d_pngPtr);

    /* unlike the example in the libpng documentation, we have *no* idea where
     * this file may have come from--so if it doesn't have a file gamma, don't
     * do any correction ("do no harm") */

	adjust_gamma(1.0);

    /* all transformations have been registered; now update info_ptr data,
     * get d_bytesPerRow and channels, and allocate image memory */

    png_read_update_info(d_pngPtr, d_infoPtr);

    d_bytesPerRow = png_get_rowbytes(d_pngPtr, d_infoPtr);   
    bands = (int)png_get_channels(d_pngPtr, d_infoPtr);
    if (bands != RGB_bands) {
    	std::cerr << "rgbImage::readPNG: (WARNING) band mismatch." << std::endl;
    }
	d_cols = d_bytesPerRow/bands;
	d_pixel.resize(d_bytesPerRow*d_rows);
	
	png_bytep* row_pointers = new png_bytep [d_rows];  
    for (i = 0; i < d_rows; i++)
        row_pointers[i] = &d_pixel[0] + i*d_bytesPerRow;

    /* now we can go ahead and just read the whole image */

    png_read_image(d_pngPtr, row_pointers);
	delete [] row_pointers;

    png_read_end(d_pngPtr, 0);
    
    fclose(fp);
}

template<typename T>
void rgbBase<T>::writePNGfile(const char* filename) {
	FILE *fp;
	const int bitsPerByte = 8;
	
	fp = fopen(filename, "wb");
	if (! fp) {
		fprintf(stderr, "RGBtoPNG: Cannot write to file %s in binary mode.\n", filename);
		exit(1);
	}
	
	/* initialize the png file structures (PNG_LIBPNG_VER_STRING is defined
	 * in png.h).
	 */
	d_pngPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      						(void*)         0,  // (void *)user_error_ptr, 
      						(png_error_ptr) 0,  // user_error_fn, 
      						(png_error_ptr) 0); // user_warning_fn);
      						
  if (! d_pngPtr) {
  	fclose(fp);
  	fprintf(stderr, "writePNGfile: Cannot create a png write structure.\n");
  	exit (1);
  }

  d_infoPtr = png_create_info_struct(d_pngPtr);
  if (! d_infoPtr) {
  	fclose(fp);
  	png_destroy_write_struct(&d_pngPtr, (png_infopp) 0);
  	fprintf(stderr, "writePNGfile: Cannot create png info structure.\n");
  	exit(1);
  }

  
  /* libpng handles errors via calls to setjmp and longjmp. Here we create
   * the user error handler with a call to setjmp.
   */ 
   if (setjmp(d_pngPtr->jmpbuf)) {
   		/* Activated only if libpng calls longjmp in the future */
   		png_destroy_write_struct(&d_pngPtr, &d_infoPtr);
   		fprintf(stderr, "writePNGfile: Error encountered in libpng.\n");
  		fclose(fp);
  		exit(1);
  }
 
  /* Tell libpng where to write the png file */
  png_init_io(d_pngPtr, fp);
  
  /* Turn the filtering off */
  png_set_filter(d_pngPtr, 0, PNG_FILTER_NONE);
  
  /* Prepare the header info */
  d_infoPtr->width         = d_cols;
  d_infoPtr->height        = d_rows;
  d_infoPtr->bit_depth     = bitsPerByte*sizeof(T);
  d_infoPtr->sig_bit.red   = d_infoPtr->bit_depth;
  d_infoPtr->sig_bit.green = d_infoPtr->bit_depth;
  d_infoPtr->sig_bit.blue  = d_infoPtr->bit_depth;
  d_infoPtr->color_type    = PNG_COLOR_TYPE_RGB;
  d_infoPtr->interlace_type = 0; /* 0 => None, 1 => Interlaced */
  
  /* Write the header */
  png_write_info(d_pngPtr, d_infoPtr);
  
  png_bytep* row_pointers = new png_bytep [d_rows];
  
  for (int i = 0; i < d_rows; i++) {
      row_pointers[i] = &d_pixel[0] + i*d_bytesPerRow;
  }
  
  /* Process the image data */
  png_write_image(d_pngPtr, row_pointers);
  delete [] row_pointers;
  
  png_write_end(d_pngPtr, d_infoPtr);
  
  /* Clean up */
  png_destroy_write_struct(&d_pngPtr, &d_infoPtr);
  d_pngPtr = 0;
  d_infoPtr = 0;
  	
fclose(fp);
}


#endif /* __PNGIMAGE_H__ */
