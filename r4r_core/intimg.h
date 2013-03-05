/*
 * intimg.h
 *
 *  Created on: Mar 30, 2012
 *      Author: jbalzer
 */

#ifndef INTIMG_H_
#define INTIMG_H_

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
	size_t Width() { return m_data.NCols(); };

	//! Image height.
	size_t Height() { return m_data.NRows(); };

	//! Computes the integral image.
	void Compute();

	//! Access to the data.
	CDenseArray<T>& Get() { return m_data; };

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
