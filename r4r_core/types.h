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

namespace R4R {

typedef CDenseVector<double> vec;
typedef CDenseVector<float> vecf;
typedef CDenseArray<double> mat;
typedef CDenseArray<float> matf;
typedef CDenseArray<unsigned char> mmat;
typedef CSparseArray<double> smat;
typedef CSparseBandedArray<double> sbmat;


}

#endif /* TYPES_H_ */
