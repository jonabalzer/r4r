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

#ifndef R4RINTIMG_H_
#define R4RINTIMG_H_

#include <assert.h>
#include <cstring>
#include <stdio.h>
#include <iostream>

#include "darray.h"
#include "rect.h"

namespace R4R {


/*! \brief integral image
 *
 * \todo Put data into a mat object with template parameter T
 *
 *
 */
template<class T>
class CIntegralImage {

public:

	//! Standard constructor.
	CIntegralImage();

	//! Standard constructor.
    CIntegralImage(size_t width, size_t height);

	//! Image width.
    size_t Width() { return m_data.NCols(); }

	//! Image height.
    size_t Height() { return m_data.NRows(); }

	//! Computes the integral image.
	void Compute();

	//! Access to the data.
    CDenseArray<T>& Get() { return m_data; }

	//! Erases the image.
	void Clear();

	/*! \brief Increases the density at a point.
	 *
	 * \detailed If the point falls between grid cells, the mass is distributed between the neighboring vertices according
	 * to the individual areas of the cell partition.
	 */
	void AddDensity(double x, double y, T val);

	/*! \brief Increases the density at a point.
	 *
	 * \detailed Non-integral locations will be rounded down.
	 *
	 */
	void AddDensityFast(double x, double y, T val);

	/*! \brief Evaluates the integral image at corners of a rectangular window around a location.
	 *
	 * \details \f$(x,y)\f$ is the center pixel. Non-integral locations will be rounded down.
	 *
	 */
	T EvaluateFast(double x, double y, double hwidth, double hheight);

	/*! \brief Evaluates the integral image at corners of a rectangular window around a location.
	 *
	 * \details \f$(x,y)\f$ is the center pixel. Non-integral locations will be interpolated bi-linearly.
	 *
	 */
	T Evaluate(double x, double y, double hwidth, double hheight);

	//! Evaluates the integral image at corners of a rectangular window around a location.
	T EvaluateFast(CRectangle<double> roi);

protected:

	CDenseArray<T> m_data;					//!< matrix holding the image data
	size_t m_computed;						//!< counts how many times CIntegralImage::Compute() has been called

	//! Maps a point that is out of bounds to its closest points on the boundary.
    void ProjectToBoundary(double& x, double& y);

	//! Checks whether a point is in bounds.
	bool InBounds(const double& x, const double& y);

	//! Gets the value of the image using bilinear interpolation.
	T Get(double x, double y);

};

}

#endif /* INTIMG_H_ */
