/*****************************************************************************/
/* grayImage.h 
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
#include "pixelArray.h"
#include <fstream.h>

typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned short ushort;

namespace PNG {
template<typename T> 
class grayBase : public pixelArray {
protected:
	std::valarray<T> d_pixel;
	png_structp d_pngPtr;
	png_infop   d_infoPtr;
	int d_pngBitDepth;
	int d_pngColorType;
public:
	grayBase(int rows, int cols) : pixelArray(rows, cols, 1) {
		assert(rows > 0);
		assert(cols > 0);
		pixelArray::d_bytesPerRow = d_cols*sizeof(T);		
		d_pixel.resize(d_rows*d_cols);		
		d_pngPtr = 0;
		d_infoPtr = 0;
	}

	grayBase() : pixelArray() {
		pixelArray::d_bands = 1;
		d_bytesPerRow = 0;
		d_pixel = 0;
		d_pngPtr = 0;
		d_infoPtr = 0;
	}
	
	~grayBase() {}
		
	grayBase(grayBase & p) {		
		d_pixel.resize(d_rows*d_cols);
		d_pixel = p.d_pixel;  // Deep copy.
		d_bytesPerRow = d_cols*sizeof(T);
		d_pngPtr = p.d_pngPtr;
		d_infoPtr = p.d_infoPtr;
	}
		
	void clear() {
		pixelArray::clear();
		d_pixel.resize(0);
		if (d_pngPtr) png_destroy_read_struct(&d_pngPtr, &d_infoPtr, 0);
	}
	
	grayBase operator=(grayBase &p) {
		d_pixel.resize(d_rows*d_cols);
		d_pixel = p.d_pixel;
		d_pngPtr = p.d_pngPtr;
		d_infoPtr = p.d_infoPtr;
	}
	
	T& graylevel(int i, int j) {
		assert(0 <= i); assert(i < d_rows);
		assert(0 <= j); assert(j < d_cols);
		return d_pixel[(i*d_cols + j)];
	}
			
    void setPixel(int i, int j, T v) {
    	assert(0 <= i); assert(i < d_rows);
		assert(0 <= j); assert(j < d_cols);
    	d_pixel[i*d_cols + j] = v;
    }
    
    void getPixel(int i, int j, T &v) {
    	v = d_pixel[i*d_cols + j];
    }
    
    T getPixel(int i, int j) {
    	if (0 <= i && i < d_rows && 0 <= j && j < d_cols) {
    		return d_pixel[i*d_cols + j];
    	} else {
    		return 0;
    	}
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
    	d_bands = 1;
    	d_pixel.resize(d_rows*d_cols);
		d_bytesPerRow = d_cols*sizeof(T);
	}
    
    double intensity() {
    	double val=0;
    	for(int i=0; i < d_rows; i++)
    		for(int j=0; j < d_cols; j++) {
    			val += static_cast<double>(d_pixel[i*d_cols + j]);
    	}
    	return val;
    }
    
    double entropy(int dr, int dc) {
    	assert(dr > 0);
    	assert(dc > 0);
    	double h=0;
    	unsigned long total = 0;
    	
    	for (int i = 0; i < d_rows; i++)
    		for (int j = 0; j < d_cols; j++) {
    			int cellIntensity = 0;
    			for (int k = i-dr/2; k < i+dr/2; k++)
    				for (int l = j-dc/2; l < j+dc/2; l++) {
    					cellIntensity += getPixel(k, l);
    			}
    			if (cellIntensity > 0) {
    				h -=  cellIntensity * std::log(static_cast<double>(cellIntensity));
    				total += cellIntensity;
    			}
    	}
    	return std::log(static_cast<double>(total)) + h/total;
   	}
    					
    		
    virtual void drawPixels() = 0; 
    	
    bool readPNGfile(const char* filename);
    void writePNGfile(const char* filename);
protected:
	void adjust_gamma(double display_exponent);
	int readpng_init(FILE *infile);
	int readpng_get_bgcolor(uchar *red, uchar *green, uchar *blue);
};

template<typename T>
class grayImage : public grayBase<T> {
	grayImage<T>() : grayBase<T>() {}
	grayImage<T>(int r, int c) : grayBase<T>(r, c) {}
};

template<>
class grayImage<GLubyte> : public grayBase<GLubyte> {
public:
	grayImage<GLubyte>() : grayBase<GLubyte>() {}
	grayImage<GLubyte>(int r, int c) : grayBase<GLubyte>(r, c) {}
	void drawPixels() {
		glDrawPixels(static_cast<GLsizei>(d_cols), 
					 static_cast<GLsizei>(d_rows), GL_LUMINANCE, GL_UNSIGNED_BYTE, &d_pixel[0]);
	}
};

template<>
class grayImage<GLshort> : public grayBase<GLshort> {
public:
	grayImage<GLshort>() : grayBase<GLshort>() {}
	grayImage<GLshort>(int r, int c) : grayBase<GLshort>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_SHORT, &d_pixel[0]);
	}
};

template<>
class grayImage<GLint> : public grayBase<GLint> {
public:
	grayImage<GLint>() : grayBase<GLint>() {}
	grayImage<GLint>(int r, int c) : grayBase<GLint>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_INT, &d_pixel[0]);
	}
};

template<>
class grayImage<GLfloat> : public grayBase<GLfloat> {
public:
	grayImage<GLfloat>() : grayBase<GLfloat>() {}
	grayImage<GLfloat>(int r, int c) : grayBase<GLfloat>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_FLOAT, &d_pixel[0]);
	}
};

template<>
class grayImage<GLdouble> : public grayBase<GLdouble> {
public:
	grayImage<GLdouble>() : grayBase<GLdouble>() {}
	grayImage<GLdouble>(int r, int c) : grayBase<GLdouble>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_DOUBLE, &d_pixel[0]);
	}
};

template<>
class grayImage<GLbyte> : public grayBase<GLbyte> {
public:
	//grayImage<GLbyte>() : grayBase<GLbyte>() {}
	//grayImage<GLbyte>(int r, int c) : grayBase<GLbyte>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_BYTE, &d_pixel[0]);
	}
};

template<>
class grayImage<GLuint> : public grayBase<GLuint> {
public:
	//grayImage<GLuint>() : grayBase<GLuint>() {}
	//grayImage<GLuint>(int r, int c) : grayBase<GLuint>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_INT, &d_pixel[0]);
	}
};

template<>
class grayImage<GLushort> : public grayBase<GLushort> {
public:
	//grayImage<GLushort>() : grayBase<GLushort>() {}
	//grayImage<GLushort>(int r, int c) : grayBase<GLushort>(r, c) {}
	void drawPixels() {
		glDrawPixels(d_cols, d_rows, GL_LUMINANCE, GL_UNSIGNED_SHORT, &d_pixel[0]);
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
void grayBase<T>::adjust_gamma(double display_exponent) {
	double gamma;
	
    if (png_get_gAMA(d_pngPtr, d_infoPtr, &gamma))
        png_set_gamma(d_pngPtr, display_exponent, gamma);
}

/* return value = 0 for success, 1 for bad sig, 2 for bad IHDR, 4 for no mem */

 bool isPNG_File(const char* path) {
   uchar sig[8];

   std::ifstream* input = new ifstream(path, ios::in | ios::binary);
   if (! input) {
	 return false;
   }

   fread(sig, 1, 8, input);
   bool status = png_check_sig(sig, 8);
   input->close();
   return status;
 }


template<typename T>
int grayBase<T>::readpng_init(FILE *infile)
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
int grayBase<T>::readpng_get_bgcolor(uchar *red, uchar *green, uchar *blue)
{
    png_color_16p pBackground;


    /* setjmp() must be called in every function that calls a PNG-reading
     * libpng function */

    if (setjmp(png_jmpbuf(d_pngPtr))) {
        png_destroy_read_struct(&d_pngPtr, &d_infoPtr, NULL);
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
bool grayBase<T>::readPNGfile(const char* filename) {
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

    if (setjmp(png_jmpbuf(d_pngPtr))) {
        png_destroy_read_struct(&d_pngPtr, &d_infoPtr, 0);
        fprintf(stderr, "readPNGfile: setjmp error");
		return false;
    }


    /* DON'T expand palette images to RGB, low-bit-depth grayscale images to 8 bits,
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
        //png_set_gray_to_rgb(d_pngPtr);

    /* unlike the example in the libpng documentation, we have *no* idea where
     * this file may have come from--so if it doesn't have a file gamma, don't
     * do any correction ("do no harm") */

	adjust_gamma(1.0);

    /* all transformations have been registered; now update info_ptr data,
     * get d_bytesPerRow and channels, and allocate image memory */

    png_read_update_info(d_pngPtr, d_infoPtr);

    d_bytesPerRow = png_get_rowbytes(d_pngPtr, d_infoPtr);   
    bands = (int)png_get_channels(d_pngPtr, d_infoPtr);
    if (bands != 1) {
    	std::cerr << "grayImage::readPNG: (WARNING) band mismatch." << std::endl;
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
	return true;
}

template<typename T>
void grayBase<T>::writePNGfile(const char* filename) {
	FILE *fp;
	const int bitsPerByte = 8;
	
	fp = fopen(filename, "wb");
	if (! fp) {
		fprintf(stderr, "grayBase<T>::writePNGfile: Cannot write to file %s in binary mode.\n", filename);
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
  	png_destroy_write_struct(&d_pngPtr, (png_infopp) NULL);
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
  d_infoPtr->sig_bit.gray   = d_infoPtr->bit_depth;
  d_infoPtr->color_type    = PNG_COLOR_TYPE_GRAY;
  d_infoPtr->interlace_type = 0; /* 0 => None, 1 => Interlaced */
  
  /* Write the header */
  png_write_info(d_pngPtr, d_infoPtr);
  
  png_bytep* row_pointers = new png_bytep [d_rows];
  for (int i = 0; i < d_rows; i++)
      row_pointers[i] = &d_pixel[0] + i*d_bytesPerRow;
        
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

} // namespace PNG

#endif /* __PNGIMAGE_H__ */
