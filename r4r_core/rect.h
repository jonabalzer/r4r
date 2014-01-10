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

#ifndef R4RRECT_H_
#define R4RRECT_H_

#include <iostream>

#include <opencv2/opencv.hpp>

#include "vecn.h"

namespace R4R {

/*! \brief rectangle
 *
 * \details The rectangle is represented as the image of the \f$[-1,1]\times[-1,1]\f$ under a similarity
 * transform. Arrays are are used to avoid possible computational overhead resulting from passing data to linear algebra
 * routines of the core library.
 *
 *
 *
 */
template<class T>
class CRectangle {

	template<class Rect, class Array> friend class CFeatureDescriptor;

public:

	//! Standard constructor creates \f$[-1,1]\times[-1,1]\f$.
	CRectangle();

	//! Parametrized constructor.
    CRectangle(T tx, T ty, T sx, T sy, T phi = 0);

	//! Transforms a point from the local coordinate system at the center of the rectangle.
    CVector<T,2> TransformFrom(T x, T y);

	//! Transforms a point from the local coordinate system at the center of the rectangle.
    CVector<T,2> TransformFrom(CVector<T,2> x) { return TransformFrom(x.Get(0),x.Get(1)); }

	//! Returns top-left corner.
    CVector<T,2> TopLeft() { return TransformFrom(-1,-1); }

	//! Returns bottom-right corner.
    CVector<T,2> BottomRight() { return TransformFrom(1,1); }

	//! Returns bottom-left corner.
    CVector<T,2> BottomLeft() { return TransformFrom(-1,1); }

	//! Returns top-right corner.
    CVector<T,2> TopRight() { return TransformFrom(1,-1); }

	//! Sets the side lengths.
    void SetSize(T hw, T hh) { m_s[0] = hw; m_s[1] = hh; }

	//! Scales the rectangle non-uniformly.
    void Scale(T sx, T sy) { m_s[0] *= sx; m_s[1] *= sy; }

	//! Scales the rectangle uniformly.
    void Scale(T s) { m_s[0] *= s; m_s[1] *= s; }

    //! Returns width. This an "oriented" width.
	T Width() { return 2*m_s[0]; }

    //! Returns height. This is an "oriented" height.
	T Height() { return 2*m_s[1]; }

	//! Computes area.
    T Area() { return 4*m_s[0]*m_s[1]; }

	//! Translates the rectangle.
    void Translate(T dx, T dy) { m_t[0] += dx; m_t[1] += dy; }

	//! Moves the rectangle to a new location.
    void MoveTo(T x, T y) { m_t[0] = x; m_t[1] = y; }

	//! Returns barycenter.
    CVector<T,2> Barycenter();

	//! Rotates the rectangle into a given orientation.
    void RotateTo(T phi) { m_phi = phi; }

	//! Rotates the rectangle by a given angle.
    void Rotate(T dphi) { m_phi += dphi; }

	//! Writes vertices of rectangle to a stream.
	template <class U> friend std::ostream& operator<<(std::ostream& os, CRectangle<U>& x);

private:

	T m_t[2];						//! translation
    T m_s[2];						//! scaling, can be negative making the transformation non-orientation preserving
    T m_phi;                        //! rotation angle

};





}

#endif /* RECT_H_ */
