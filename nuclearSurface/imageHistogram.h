#ifndef _IMAGE_HISTOGRAM_H_
#define _IMAGE_HISTOGRAM_H_

//
//  imageHistogram.h
//  lsmModeler
//
//  Created by Robert Snapp on 2007-02-06.
//  Copyright (c) 2007 Robert R. Snapp. All rights reserved.
//

#include "env.h"  // for OpenGL headers.
#include "grayImage3d.h"
#include "lsmVoxelImage.h"
#include <algorithm>
#include <valarray>
#include <cmath>

class ImageHistogram {
private:
	std::valarray<ulong> d_bin;
	double d_mean; // mean
	double d_sdev; // standard deviation
	
public:
	ImageHistogram(ulong nBins = 1) {
		nBins = std::max<ulong>(1, nBins);
		d_bin.resize(nBins, 0);
	}
	
	ImageHistogram(const ImageHistogram &hist) {
		d_bin.resize(hist.d_bin.size());
		for(ulong i = 0; i < hist.d_bin.size(); i++) {
			d_bin[i] = hist.d_bin[i];
		}
		d_mean = hist.d_mean;
		d_sdev = hist.d_sdev;
	}
	
	ImageHistogram& operator=(const ImageHistogram& hist) {
		d_bin.resize(hist.d_bin.size());
		for(ulong i = 0; i < hist.d_bin.size(); i++) {
			d_bin[i] = hist.d_bin[i];
		}
		d_mean = hist.d_mean;
		d_sdev = hist.d_sdev;
		return *this;
	}
	
	void processImage(const grayImage3d<GLubyte>& image) {
		double delta = 256.0/d_bin.size();
		ulong sum = 0;
		ulong count = 0;
		ulong sum2 = 0;

		for (ulong l = 0; l < image.layers(); l++ )
			for (ulong i = 0; i < image.rows(); i++)
				for (ulong j = 0; j < image.cols(); j++) {
					ulong val = image.getVoxelClip(i, j, l);
					ulong binIndex = val/delta;
					d_bin[binIndex]++;
					sum  += val;
					sum2 += val*val;
					count++;
		}

		d_mean = static_cast<double>(sum)/count;
		d_sdev = std::sqrt(static_cast<double>(sum2)/count - d_mean*d_mean);
		return;
	}
	
	void processImage(const lsmVoxelImage& image, ulong channel) {
		double delta = 256.0/d_bin.size();
		ulong sum = 0;
		ulong count = 0;
		ulong sum2 = 0;

		for (ulong l = 0; l < image.layers(); l++ )
			for (ulong i = 0; i < image.rows(); i++)
				for (ulong j = 0; j < image.cols(); j++) {
					ulong val = image.getVoxelClip(i, j, l, channel);
					ulong binIndex = val/delta;
					d_bin[binIndex]++;
					sum  += val;
					sum2 += val*val;
					count++;
		}

		d_mean = static_cast<double>(sum)/count;
		d_sdev = std::sqrt(static_cast<double>(sum2)/count - d_mean*d_mean);
		return;
	}
	
	void display(bool autoScale = true, float sx = 1.0, float sy = 1.0) const {
		glPushMatrix();
		if (autoScale) {
		  sx = 1.0/d_bin.size();
		  sy = 1.0/d_bin.max();
		}
		glScalef(sx, sy, 1.0);	
		glBegin(GL_QUADS);
			for (ulong i = 0; i < d_bin.size(); i++) {
				glVertex2i(i,   0);
				glVertex2i(i+1, 0);
				glVertex2i(i+1, static_cast<GLint>(d_bin[i]));
				glVertex2i(i,   static_cast<GLint>(d_bin[i]));
			}
		glEnd();
		glPopMatrix();	
	}

	double getMean() const {return d_mean;}
	double getSDev() const {return d_sdev;}
	ulong size() const {return d_bin.size();}
	ulong operator[](ulong i) const { return d_bin[i];}
	
	void clear() {
		d_bin = 0;
		d_mean = 0;
		d_sdev = 0;
	}

};

#endif /* _IMAGE_HISTOGRAM_H_ */
