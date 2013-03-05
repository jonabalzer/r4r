/*
 * interp.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: jbalzer
 */



#include "interp.h"

using namespace std;
using namespace cv;

namespace R4R {

vector<Mat> CImageInterpolation::BuildGaussianPyramidForDetection(const Mat& img, size_t nlevels, size_t filtersize, double sigma) {

	vector<Mat> pyramid;

	Mat p0 = img.clone();

	GaussianBlur(img,p0,Size(filtersize,filtersize),sigma);

	pyramid.push_back(p0);

	for(size_t i=0; i<nlevels; i++) {

		Mat dst = pyramid.back().clone();

		GaussianBlur(pyramid.back(),dst,Size(filtersize,filtersize),sigma);

		Mat dst2 = CImageInterpolation::Downsample(dst);

		pyramid.push_back(dst2);

	}

	return pyramid;

}



}



