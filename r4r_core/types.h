/*
 * types.h
 *
 *  Created on: May 11, 2012
 *      Author: jbalzer
 */

#ifndef TYPES_H_
#define TYPES_H_

#include "darray.h"
#include "sarray.h"
#include <complex>

namespace R4R {


typedef CDenseVector<double> vec;
typedef CDenseVector<std::complex<double> > vecc;
typedef CDenseVector<float> vecf;
typedef CDenseArray<double> mat;
typedef CDenseArray<float> matf;
typedef CDenseArray<unsigned char> mmat;
typedef CSparseArray<double> smat;
typedef CSparseBandedArray<double> sbmat;
typedef CNVector<double,3> vec3;
typedef CNVector<unsigned char,3> rgb;


}

#endif /* TYPES_H_ */
