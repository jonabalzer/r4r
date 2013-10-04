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

#ifndef R4RTYPES_H_
#define R4RTYPES_H_

#include "vecn.h"
#include "darray.h"
#include "sarray.h"
#include <complex>

namespace R4R {

typedef CVector<double,3> vec3;
typedef CVector<float,3> vec3f;
typedef CVector<double,2> vec2;
typedef CVector<float,2> vec2f;
typedef CVector<unsigned char,3> rgb;


typedef CDenseVector<double> vec;
typedef CDenseVector<std::complex<double> > vecc;
typedef CDenseVector<float> vecf;


typedef CDenseArray<double> mat;
typedef CDenseArray<float> matf;
typedef CDenseArray<unsigned char> mmat;

typedef CSparseArray<double> smat;
typedef CSparseArray<float> smatf;
typedef CSparseBandedArray<double> sbmat;

}

#endif /* TYPES_H_ */
