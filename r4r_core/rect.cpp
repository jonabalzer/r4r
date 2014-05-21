//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
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
//////////////////////////////////////////////////////////////////////////////////

#include "rect.h"
#include <iostream>
#include <assert.h>
#include <math.h>

#define MARGIN 2

using namespace std;

namespace R4R {

template <class T>
CRectangle<T>::CRectangle() {

	m_t[0] = 0;
	m_t[1] = 0;
	m_s[0] = 1;
	m_s[1] = 1;
	m_phi = 0;

}


template <class T>
CRectangle<T>::CRectangle(T tx, T ty, T sx, T sy, T phi) {

	m_t[0] = tx;
	m_t[1] = ty;
	m_s[0] = sx;
	m_s[1] = sy;
	m_phi = phi;

}

template <class T>
CVector<T,2> CRectangle<T>::TransformFrom(T x, T y) {

    CVector<T,2> result;

	if(m_phi!=0) {

		T cphi = (T)cos(m_phi);
		T sphi = (T)sin(m_phi);

		result(0) = m_s[0]*(cphi*x - sphi*y) + m_t[0];
		result(1) = m_s[1]*(sphi*x + cphi*y) + m_t[1];

	} else {

		result(0) = m_s[0]*x + m_t[0];
		result(1) = m_s[1]*y + m_t[1];

	}

	return result;

}


template <class T>
CVector<T,2> CRectangle<T>::Barycenter() {

    return  { m_t[0], m_t[1] };

}

template <class U>
ostream& operator<<(ostream& os, CRectangle<U>& x) {

    os << "[ ";
    os << x.TopLeft() << endl;
    os << x.BottomLeft() << endl;
    os << x.BottomRight() << endl;
    os << x.TopRight() << endl;
    os << " ]";

	return os;

}

template class CRectangle<float>;
template class CRectangle<double>;
template class CRectangle<int>;              // this does not make sense for unsigned values!!!
template ostream& operator<< (ostream& os, CRectangle<double>& x);
template ostream& operator<< (ostream& os, CRectangle<float>& x);

}
