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

#include "intimg.h"
#include <math.h>
#include <iostream>

using namespace std;

namespace R4R {

template <class T>
CIntegralImage<T>::CIntegralImage():
	m_data(0,0),
	m_computed(0) {


}

template <class T>
CIntegralImage<T>::CIntegralImage(size_t width, size_t height):
	m_data(height,width),
	m_computed(0) {

}

template <class T>
void CIntegralImage<T>::Clear() {

	m_data.Scale(0);

}

template <class T>
void CIntegralImage<T>::Compute() {

	if(m_data.NCols()>=m_data.NRows()) {

		for(size_t d=1; d<m_data.NRows(); d++) {

			// process the quadratic part of the matrix
			for(size_t i=0; i<d; i++) {

				if(i==0) {

					// column above (d+d)
					m_data(d,i) += m_data(d-1,i);

					// row left of (d+d)
					m_data(i,d) += m_data(i,d-1);

				}
				else {

					// column above (d+d)
					m_data(d,i) += m_data(d,i-1) + m_data(d-1,i) - m_data(d-1,i-1);

					// row left of (d+d)
					m_data(i,d) += m_data(i,d-1) + m_data(i-1,d) - m_data(i-1,d-1);

				}

			}

			// value at (d+d)
			m_data(d,d) += m_data(d,d-1) + m_data(d-1,d) - m_data(d-1,d-1);

		}

		// process the remaining columns
		for(size_t j=m_data.NRows(); j<m_data.NCols(); j++) {

			for(size_t i=0; i<m_data.NRows(); i++) {

				if(i==0)
					m_data(i,j) += m_data(i,j-1);
				else
					m_data(i,j) += m_data(i,j-1) + m_data(i-1,j) - m_data(i-1,j-1);

			}

		}

	}
	else {

		for(size_t d=1; d<m_data.NCols(); d++) {

			for(size_t i=0; i<d; i++) {

				if(i==0) {

					// column above (d+d)
					m_data(d,i) += m_data(d-1,i);

					// row left of (d+d)
					m_data(i,d) += m_data(i,d-1);

				}
				else {

					// column above (d+d)
					m_data(d,i) += m_data(d,i-1) + m_data(d-1,i) - m_data(d-1,i-1);

					// row left of (d+d)
					m_data(i,d) += m_data(i,d-1) + m_data(i-1,d) - m_data(i-1,d-1);

				}


			}

			// value at (d+d)
			m_data(d,d) += m_data(d,d-1) + m_data(d-1,d) - m_data(d-1,d-1);

		}

		// process the remaining row
		for(size_t i=m_data.NCols(); i<m_data.NRows(); i++) {

			for(size_t j=0; j<m_data.NCols(); j++) {

				if(j==0)
					m_data(i,j) += m_data(i-1,j);
				else
					m_data(i,j) += m_data(i,j-1) + m_data(i-1,j) - m_data(i-1,j-1);


			}

		}

	}

	m_computed++;

}

template<>
void CIntegralImage<double>::AddDensity(double x, double y, double val) {

	if(InBounds(x,y)) {

		double dx, dy;
		dx = x - floor(x);
		dy = y - floor(y);

		size_t i, j;
		i = size_t(y);
		j = size_t(x);

		m_data(i,j) += dx*dy;
		m_data(i+1,j) += dx*(1-dy);
		m_data(i,j+1) += dy*(1-dx);
		m_data(i+1,j+1) += (1-dx)*(1-dy);

	}

}

template <class T>
void CIntegralImage<T>::AddDensityFast(double x, double y, T val) {

	if(InBounds(x,y)) {

		size_t i, j;
		i = size_t(y);
		j = size_t(x);

		m_data(i,j) += val;

	}

}



template <class T>
void CIntegralImage<T>::ProjectToBoundary(double& x, double& y) {

	if(x<0)
		x = 0;

	if(x>=m_data.NCols())
		x = m_data.NCols() - 1;

	if(y<0)
		y = 0;

	if(y>=m_data.NRows())
		y = m_data.NRows() - 1;

}

template <class T>
bool CIntegralImage<T>::InBounds(const double& x, const double& y) {

	return x>=0 && y>=0 && x<m_data.NCols() && y<m_data.NRows();

}

template <class T>
T CIntegralImage<T>::EvaluateFast(double x, double y, double hwidth, double hheight) {

	assert(m_computed>0);

	double xl, xu, yl, yu;

	xl = x - hwidth;
	xu = x + hwidth;
	yl = y - hheight;
	yu = y + hheight;

	ProjectToBoundary(xl,yl);
	ProjectToBoundary(xu,yu);

	size_t il, jl, iu, ju;
	il = (size_t)yl;
	jl = (size_t)xl;
	iu = (size_t)yu;
	ju = (size_t)xu;

	return m_data(il,jl) + m_data(iu,ju) - m_data(iu,jl) - m_data(il,ju);

}

template <class T>
T CIntegralImage<T>::Evaluate(double x, double y, double hwidth, double hheight) {

	assert(m_computed>0);

	double xl, xu, yl, yu;

	xl = x - hwidth;
	xu = x + hwidth;
	yl = y - hheight;
	yu = y + hheight;

	ProjectToBoundary(xl,yl);
	ProjectToBoundary(xu,yu);

	return Get(xl,yl) + Get(xu,yu) - Get(xu,yl) - Get(xl,yu);

}

template <class T>
T CIntegralImage<T>::EvaluateFast(CRectangle<double> roi) {

	assert(m_computed>0);

	vec tl = roi.TopLeft();
	vec br = roi.BottomRight();

	ProjectToBoundary(tl(0),tl(1));
	ProjectToBoundary(br(0),br(1));

	size_t il, jl, iu, ju;
	il = (size_t)tl(1);
	jl = (size_t)tl(0);
	iu = (size_t)br(1);
	ju = (size_t)br(0);

	return m_data(il,jl) + m_data(iu,ju) - m_data(iu,jl) - m_data(il,ju);


}

template <class T>
T CIntegralImage<T>::Get(double x, double y) {

	// FIXME: implement interpolation

	return 0;


}

template class CIntegralImage<size_t>;



}
