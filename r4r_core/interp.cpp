/*////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the R4R library.
//
// The R4R library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The R4R library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////*/

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



