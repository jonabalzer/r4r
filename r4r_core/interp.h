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

#ifndef R4RINTERP_H_
#define R4RINTERP_H_

#include "opencv2/opencv.hpp"
#include "types.h"

using namespace std;

namespace R4R {

/*! \brief image interpolation
 *
 *
 *
 */
class CImageInterpolation {

public:

	//! Interpolates an image bilinearly.
	template <class T = unsigned char>
	static inline double Bilinear(const cv::Mat& img, const double u, const double v) {

		int i = (int)floor(v);
		int j = (int)floor(u);

		if(i<0 || i>=img.rows || j<0 || j>=img.cols)
			return 0;

		double vd = v - i;
		double ud = u - j;

		double I00, I01, I10, I11, I0, I1;
		I00 = img.at<T>(i,j);
		I01 = img.at<T>(i,j+1);
		I10 = img.at<T>(i+1,j);
		I11 = img.at<T>(i+1,j+1);
		I0 = (1-ud)*I00 + ud*I01;
		I1 = (1-ud)*I10 + ud*I11;

		return (1-vd)*I0 + vd*I1;

	}

	//! Computes the gradient of an image using bilinear interpolation.
	template <class T = unsigned char>
	static inline double Gradient(const cv::Mat& img, const double u, const  double v, const bool dir) {

		double I1, I0;

		if(dir) {


			I1 = Bilinear<T>(img,u+1.0,v);
			I0 = Bilinear<T>(img,u-1.0,v);


		} else {


			I1 = Bilinear<T>(img,u,v+1.0);
			I0 = Bilinear<T>(img,u,v-1.0);

		}

		return 0.5*(I1-I0);

	}


	//! Interpolates a matrix bilinearly.
	static inline double Bilinear(mat& A, double u, double v) {

		if(u<0 || u>=(double)(A.NRows()-1) || v<0 || v>=(double)(A.NCols()-1))
			return 0;

		int i = (int)floor(v);
		int j = (int)floor(u);

		double vd = v - i;
		double ud = u - j;

		double A00, A01, A10, A11, A0, A1;
		A00 = A.Get(i,j);
		A01 = A.Get(i,j+1);
		A10 = A.Get(i+1,j);
		A11 = A.Get(i+1,j+1);
		A0 = (1-ud)*A00 + ud*A01;
		A1 = (1-ud)*A10 + ud*A11;

		return (1-vd)*A0 + vd*A1;

	}

	template <class T = size_t>
	static inline void ProjectRect(T& utl, T& vtl, T& w, T& h, size_t width, size_t height) {

		// project top-left corner
		if(utl<0)
			utl=0;

		if(utl>=width)
			utl = width-1;

		if(vtl<0)
			vtl = 0;

		if(vtl>=height)
			vtl = height - 1;

		if(utl+w>width)			// br corner may be ouside the image
			w = width - utl;

		if(vtl+h>height)
			h = height - vtl;

	}

	static std::vector<cv::Mat> BuildGaussianPyramidForDetection(const cv::Mat& img, size_t nlevels, size_t filtersize, double sigma);

	template <class T = unsigned char>
	static inline cv::Mat Downsample(const cv::Mat& img) {

		size_t nrows = img.rows/2;
		size_t ncols = img.cols/2;

		cv::Mat out(cv::Size(ncols,nrows),img.type());

		for(size_t i=0; i<nrows; i++) {

			for(size_t j=0; j<ncols; j++) {

				out.at<T>(i,j) = img.at<T>(2*i,2*j);

			}

		}

		return out;

	}


};





}

#endif /* INTERP_H_ */
