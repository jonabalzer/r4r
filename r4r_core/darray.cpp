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

#include "darray.h"
#include "types.h"

#include <string.h>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <chrono>
#include <random>


#ifdef _OPENMP
#include <omp.h>
#include <parallel/algorithm>
#endif

using namespace std;

namespace R4R {

template <class T>
CDenseArray<T>::CDenseArray():
	m_nrows(0),
	m_ncols(0),
	m_transpose(false),
    m_data() {

}

template <class T>
CDenseArray<T>::CDenseArray(size_t nrows, size_t ncols, T val):
	m_nrows(nrows),
	m_ncols(ncols),
    m_transpose(false),
#ifndef __SSE4_1__
    m_data(new T[nrows*ncols]) {
#else
    m_data((T*)_mm_malloc(nrows*ncols*sizeof(T),16),CDenseMatrixDeallocator<T>()){
#endif

    fill_n(m_data.get(),nrows*ncols,val);

}

template <class T>
CDenseArray<T>::CDenseArray(const CDenseArray& array):
	m_nrows(array.m_nrows),
	m_ncols(array.m_ncols),
    m_transpose(array.m_transpose),
    m_data(array.m_data) {

}

template <class T>
CDenseArray<T> CDenseArray<T>::operator=(const CDenseArray<T>& array) {

    if(this==&array)
        return *this;

    m_nrows = array.m_nrows;
    m_ncols = array.m_ncols;
    m_transpose = array.m_transpose;
    m_data = array.m_data;

    return *this;

}

template <class T>
CDenseArray<T>::CDenseArray(size_t nrows, size_t ncols, shared_ptr<T> data):
	m_nrows(nrows),
	m_ncols(ncols),
    m_transpose(false),
    m_data(data) {}

template <class T>
void CDenseArray<T>::Set(shared_ptr<T> data) {

    m_data = data;

}

template <class T>
CDenseArray<T> CDenseArray<T>::Clone() {

    CDenseArray<T> result(m_nrows,m_ncols);

    // copy data
    T* newdata = result.Data().get();
    T* thisdata = m_data.get();
    memcpy(newdata,thisdata,m_nrows*m_ncols*sizeof(T));

    return result;

}

template <class T>
CDenseArray<T>::~CDenseArray() {

    // nothing to do with smart pointers

}

template <class T>
void CDenseArray<T>::Eye() {

	this->Scale(0);

	for(size_t i=0; i<min(this->NRows(),this->NCols());i++)
		this->operator ()(i,i) = 1;

}


template <>
void CDenseArray<float>::Rand(float min, float max) {

    unsigned short seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    uniform_real_distribution<float> distribution(min,max);

    for(size_t i=0; i<m_nrows;i++) {

        for(size_t j=0; j<m_ncols; j++) {

            this->operator ()(i,j) = distribution(generator);

        }

    }

}

template <>
void CDenseArray<double>::Rand(double min, double max) {

    unsigned short seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    uniform_real_distribution<double> distribution(min,max);

    for(size_t i=0; i<m_nrows;i++) {

        for(size_t j=0; j<m_ncols; j++) {

            this->operator ()(i,j) = distribution(generator);

        }

    }

}

template <>
void CDenseArray<int>::Rand(int min, int max) {

    unsigned short seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    uniform_int_distribution<int> distribution(min,max);

    for(size_t i=0; i<m_nrows;i++) {

        for(size_t j=0; j<m_ncols; j++) {

            this->operator ()(i,j) = distribution(generator);

        }

    }

}

template <>
void CDenseArray<float>::RandN(float mu, float sigma) {

    unsigned short seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<float> distribution(mu,sigma);

    for(size_t i=0; i<m_nrows;i++) {

        for(size_t j=0; j<m_ncols; j++) {

            this->operator ()(i,j) = distribution(generator);

        }

    }

}

template <>
void CDenseArray<double>::RandN(double mu, double sigma) {

    unsigned short seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(mu,sigma);

    for(size_t i=0; i<m_nrows;i++) {

        for(size_t j=0; j<m_ncols; j++) {

            this->operator ()(i,j) = distribution(generator);

        }

    }

}

template <class T>
void CDenseArray<T>::Ones() {

    T* pdata = m_data.get();

    fill_n(pdata,NElems(),1);

}

template <class U>
ostream& operator<< (ostream& os, const CDenseArray<U>& x) {

	os.setf(ios::scientific,ios::floatfield);
	os.precision(5);

	for(size_t i=0; i<x.m_nrows; i++) {

		for(size_t j=0; j<x.m_ncols-1; j++) {

			os << (float)x.Get(i,j) << " ";

		}

		if(i<x.m_nrows-1)
			os << (float)x.Get(i,x.m_ncols-1) << endl;
		else
			os << (float)x.Get(i,x.m_ncols-1);

	}

	return os;

}

ofstream& operator<< (ofstream& os, const CDenseArray<bool>& x) {

    os << x.NRows() << " " << x.NCols() << " B1U" << endl;

    os.write((char*)(x.m_data.get()),sizeof(bool)*x.NElems());

    return os;

}

ofstream& operator<< (ofstream& os, const CDenseArray<int>& x) {

    os << x.NRows() << " " << x.NCols() << " I4S" << endl;

    os.write((char*)(x.m_data.get()),sizeof(int)*x.NElems());

    return os;

}

ofstream& operator<< (ofstream& os, const CDenseArray<float>& x) {

    os << x.NRows() << " " << x.NCols() << " F4S" << endl;

    os.write((char*)(x.m_data.get()),sizeof(float)*x.NElems());

    return os;

}

ofstream& operator<< (ofstream& os, const CDenseArray<size_t>& x) {

    os << x.NRows() << " " << x.NCols() << " L8U" << endl;

    os.write((char*)(x.m_data.get()),sizeof(size_t)*x.NElems());

    return os;

}

ofstream& operator<< (ofstream& os, const CDenseArray<double>& x) {

    os << x.NRows() << " " << x.NCols() << " D8S" << endl;

    os.write((char*)(x.m_data.get()),sizeof(double)*x.NElems());

    return os;

}

template<>
ETYPE CDenseArray<bool>::GetType() { return ETYPE::B1U; }

template<>
ETYPE CDenseArray<int>::GetType() { return ETYPE::I4S; }

template<>
ETYPE CDenseArray<size_t>::GetType() { return ETYPE::L8U; }

template<>
ETYPE CDenseArray<float>::GetType() { return ETYPE::F4S; }

template<>
ETYPE CDenseArray<double>::GetType() { return ETYPE::D8S; }

template<class U>
istream& operator >> (istream& is, CDenseArray<U>& x) {

	for(size_t i=0; i<x.m_nrows; i++) {

		for(size_t j=0; j<x.m_ncols; j++) {

			 is >> x(i,j);

		}

	}

	return is;

}

template <class U>
ifstream& operator >> (ifstream& in, CDenseArray<U>& x) {

    // read information
    size_t nrows, ncols;
    in >> nrows;
    in >> ncols;

    // resize storage if necessary
    if(nrows!=x.NRows() || ncols!=x.NCols()) {

        x.m_data.reset();

#ifndef __SSE4_1__
        x.m_data.reset(new U[nrows*ncols]);
#else
        x.m_data.reset((U*)_mm_malloc(nrows*ncols*sizeof(U),16));
#endif
        x.m_nrows = nrows;
        x.m_ncols = ncols;
        x.m_transpose = false;

    }

    char sign, type;
    size_t nbytes;

    in >> type;
    in >> nbytes;
    in >> sign;
    in.get();

    stringstream ss;
    ss << type << nbytes << sign;

    // create buffer
    char* buffer = new char[nrows*ncols*nbytes];

    // read data
    in.read(buffer,nbytes*nrows*ncols);

    U* pdata = x.m_data.get();

    if(ss.str()=="B1U") {

        bool* cbuffer = (bool*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="S2S") {

        short* cbuffer = (short*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="S2U") {

        unsigned short* cbuffer = (unsigned short*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="C4S") {

        char* cbuffer = (char*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="C4U") {

        unsigned char* cbuffer = (unsigned char*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="I4S") {

        int* cbuffer = (int*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="I4U") {

        unsigned int* cbuffer = (unsigned int*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="F4S") {

        float* cbuffer = (float*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }


    if(ss.str()=="L8S") {

        long* cbuffer = (long*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="L8U") {

        unsigned long* cbuffer = (unsigned long*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    if(ss.str()=="D8S") {

        double* cbuffer = (double*)buffer;

        // copy data and cast
        for(size_t i=0; i<nrows*ncols; i++)
            pdata[i] = (U)cbuffer[i];

    }

    delete [] buffer;

    return in;

}

template <class T>
bool CDenseArray<T>::Normalize() {

    double norm = Norm2();

	if(norm>0) {

        T s = T(1/norm);

		Scale(s);

		return 0;

	}

	return 1;

}


template <class T>
T CDenseArray<T>::Trace() const {

	size_t m = min(m_nrows,m_ncols);

	T sum = 0;

	for(size_t i=0; i<m; i++)
		sum += Get(i,i);

	return sum;

}

template <class T>
T CDenseArray<T>::Determinant() const {

	assert((m_ncols==2 && m_nrows==2) || (m_ncols==3 && m_nrows==3));

	T det = 0;

	switch (m_nrows) {

		case 2:

			det = Get(0,0)*Get(1,1) - Get(0,1)*Get(1,0);

			break;

		case 3:

			det = Get(0,0)*(Get(1,1)*Get(2,2)-Get(1,2)*Get(2,1)) - Get(0,1)*(Get(1,0)*Get(2,2)-Get(2,0)*Get(1,2)) + Get(0,2)*(Get(1,0)*Get(2,1)-Get(2,0)*Get(1,1));

			break;

		default:

			det = 0;

			break;
	}

	return det;

}


template <class T>
void CDenseArray<T>::Transpose() {

	size_t temp = m_nrows;

	m_nrows = m_ncols;
	m_ncols = temp;

	m_transpose = !m_transpose;

}

template <class T>
double CDenseArray<T>::Norm2() const {

    CMercerKernel<T> kernel(m_nrows*m_ncols);

    T* pdata = m_data.get();

    return kernel.ComputeKernelNorm(pdata);

}

template<>
double CDenseArray<bool>::HammingNorm() {

    double sum = 0;

    bool* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        sum += fabs((double)pdata[i]);

    return sum;

}

template <class T>
double CDenseArray<T>::Norm1() const {

    CHellingerKernel<T> kernel(m_nrows*m_ncols);

    T* pdata = m_data.get();

    return kernel.ComputeKernelNorm(pdata);

}

template <class T>
double CDenseArray<T>::Norm(double p) const {

    assert(p>0);

	double sum = 0;

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        sum += pow(fabs((double)pdata[i]),p);

    return pow(sum,1/p);

}

template <class T>
T CDenseArray<T>::Sum() const {

	T sum = 0;

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        sum += pdata[i];

	return sum;

}


template <class T>
void CDenseArray<T>::Abs() {

    T* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = fabs(pdata[i]);

}


template <class T>
T CDenseArray<T>::Get(size_t i, size_t j) const {

    assert(i<m_nrows && j<m_ncols);

    T* pdata = m_data.get();

	if(!m_transpose)
        return pdata[m_nrows*j + i];
	else
        return pdata[m_ncols*i + j];

}

template <class T>
CDenseVector<T> CDenseArray<T>::GetColumn(size_t j) const {

    if(m_transpose)
        return GetRow(j);

    // create shared pointer to col data
    shared_ptr<T> colptr(m_data,m_data.get()+j*m_nrows);

    // create column view
    CDenseVector<T> col(m_nrows,colptr);

	return col;

}

template <class T>
void CDenseArray<T>::SetColumn(size_t j, const CDenseVector<T>& col) {

    assert(col.NCols()==1 && col.NRows() == NRows() && j<=NCols());

	for(size_t i=0; i<col.NRows(); i++)
		this->operator ()(i,j) = col.Get(i);

}

template <class T>
void CDenseArray<T>::SetRow(size_t i, const CDenseVector<T>& row) {

    assert(row.NRows()==1 && row.NCols() == NCols() && i<=NRows());

    for(size_t j=0; j<row.NCols(); j++)
        this->operator ()(i,j) = row.Get(j);

}

template <class T>
CDenseVector<T> CDenseArray<T>::GetRow(size_t i) const {

    if(m_transpose)
        return GetColumn(i);

    // since matrix is stored in col-major order, we have to make a hard-copy
    CDenseVector<T> row(m_ncols);
    row.Transpose();

	for(size_t j=0; j<m_ncols; j++)
		row(j) = Get(i,j);

	return row;

}

template <class T>
T& CDenseArray<T>::operator()(size_t i, size_t j) {

    assert(i<m_nrows && j<m_ncols);

    T* pdata = m_data.get();

	if(!m_transpose)
        return pdata[m_nrows*j + i];
	else
        return pdata[m_nrows*i + j];				// m_nrows is replaced by m_ncols in transpose method

}

template <class T>
CDenseArray<T> CDenseArray<T>::operator+(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = pdata[i] + scalar;

	return result;

}

template <class T>
CDenseArray<T> CDenseArray<T>::operator-(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();
    T* pdatares = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = pdata[i] - scalar;

	return result;

}


template <class T>
CDenseArray<T> CDenseArray<T>::operator+(const CDenseArray& array) const {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

    CDenseArray<T> result(array.m_nrows,array.m_ncols);

    // this does not work low-level because of possible transpositions
	for(size_t i=0; i<m_nrows; i++) {

		for(size_t j=0; j<m_ncols; j++) {

			result(i,j) = this->Get(i,j) + array.Get(i,j);

		}

	}

	return result;

}


template <class T>
CDenseArray<T> CDenseArray<T>::operator^(const CDenseArray& array) const {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

    CDenseArray<T> result(array.m_nrows,array.m_ncols);

    T* px = m_data.get();
    T* py = array.m_data.get();
    T* pz = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pz[i] = px[i]*py[i];

	return result;

}

template <class T>
double CDenseArray<T>::InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y) {

	assert(x.m_nrows==y.m_nrows && x.m_ncols==y.m_ncols);

    T* px = x.m_data.get();
    T* py = y.m_data.get();

    CMercerKernel<T> kernel(x.m_nrows*x.m_ncols);

    return kernel.Evaluate(px,py);

}

template <class T>
double CDenseArray<T>::InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y, CMercerKernel<T>& kernel) {

    assert(x.m_nrows==y.m_nrows && x.m_ncols==y.m_ncols);

    T* px = x.m_data.get();
    T* py = y.m_data.get();

    return kernel.Evaluate(px,py);

}

template <class T>
CDenseArray<T> CDenseArray<T>::KroneckerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y) {

    CDenseArray<T> result(x.NRows()*y.NRows(),x.NCols()*y.NCols());

	for(size_t i=0; i<x.NRows(); i++) {

		for(size_t j=0; j<x.NCols(); j++) {

			for(size_t k=0; k<y.NRows(); k++) {

				for(size_t l=0; l<y.NCols(); l++) {

					result(i*y.NRows()+k,j*y.NCols()+l) = x.Get(i,j)*y.Get(k,l);

				}

			}

		}

	}

	return result;

}

template <class T>
CDenseArray<T> CDenseArray<T>::Transpose(const CDenseArray<T>& x) {

	CDenseArray<T> result(x);

	result.Transpose();

	return result;

}


template <class T>
CDenseArray<T> CDenseArray<T>::operator-(const CDenseArray<T>& array) const {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

    CDenseArray<T> result(array.m_nrows,array.m_ncols);

	for(size_t i=0; i<m_nrows; i++) {

		for(size_t j=0; j<m_ncols; j++) {

			result(i,j) = this->Get(i,j) - array.Get(i,j);

		}

	}

	return result;

}

template <class T>
CDenseArray<T> CDenseArray<T>::operator*(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdatares = result.m_data.get();
    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = scalar*pdata[i];

	return result;

}

template <class T>
CDenseArray<T> CDenseArray<T>::operator/(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdatares = result.m_data.get();
    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = pdata[i]/scalar;

	return result;

}

template <class T>
CDenseArray<T> CDenseArray<T>::operator*(const CDenseArray<T>& array) const {

	assert(m_ncols==array.m_nrows);

    CDenseArray<T> result(m_nrows,array.m_ncols);

#pragma omp parallel for
	for(size_t i=0; i<m_nrows; i++) {

		for(size_t j=0; j<array.m_ncols; j++) {

			T sum = 0;

			for(size_t k=0; k<m_ncols; k++)
				sum += Get(i,k)*(array.Get(k,j));

			result(i,j) = sum;

		}

	}


	return result;

}

template <class T>
CDenseVector<T> CDenseArray<T>::operator*(const CDenseVector<T>& vector) const {

	assert(m_ncols==vector.m_nrows);

    CDenseVector<T> result(m_nrows);

#pragma omp parallel for
	for(size_t i=0; i<m_nrows; i++) {

		T sum = 0;

		for(size_t k=0; k<m_ncols; k++)
			sum += Get(i,k)*(vector.Get(k));

		result(i) = sum;

	}

	return result;

}



template <class T>
void CDenseArray<T>::Scale(T scalar) {

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] *= scalar;

}

template <class T>
void CDenseArray<T>::Add(const T& scalar) {

    T* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] += scalar;

}

template <class T>
T CDenseArray<T>::Median() {

    CDenseArray<T> temp = this->Clone();

#ifdef _OPENMP
    __gnu_parallel::sort(temp.m_data.get(),temp.m_data.get()+temp.NElems());
#else
    sort(temp.m_data.get(),temp.m_data.get()+temp.NElems());
#endif

	if(temp.NElems()%2==1)
        return temp.m_data.get()[(size_t)((temp.NElems()+1)/2)-1];
	else
        return 0.5*(temp.m_data.get()[temp.NElems()/2-1]+temp.m_data.get()[temp.NElems()/2]);

}

template <class T>
T CDenseArray<T>::Variance() {

	T mean = this->Mean();

    CDenseArray<T> temp = this->Clone();

    T* pdata = temp.m_data.get();
    T* pthisdata = m_data.get();


	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = (pthisdata[i] - mean)*(pthisdata[i] - mean);

	return temp.Mean();

}

template <class T>
T CDenseArray<T>::MAD() {

    T median = Median();

    CDenseArray<T> temp = this->Clone();

    T* pdata = temp.m_data.get();
    T* pthisdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = fabs(pthisdata[i]-median);

    return temp.Median();

}

template <class T>
T CDenseArray<T>::Min() const {

	T min = numeric_limits<T>::max();

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++) {

        if(pdata[i]<=min)
            min = pdata[i];

	}

	return min;

}

template <class T>
T CDenseArray<T>::Max() const {

	T max = numeric_limits<T>::min();

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++) {

        if(pdata[i]>=max)
            max = pdata[i];

	}

	return max;

}


template class CDenseArray<float>;
template class CDenseArray<double>;
template class CDenseArray<int>;
template class CDenseArray<size_t>;
template class CDenseArray<bool>;
//template class CDenseArray<rgb>;

template ostream& operator<< (ostream& os, const CDenseArray<double>& x);
template istream& operator>> (istream& is, CDenseArray<double>& x);
template ifstream& operator>> (ifstream& is, CDenseArray<double>& x);
template ostream& operator<< (ostream& os, const CDenseArray<float>& x);
template istream& operator>> (istream& is, CDenseArray<float>& x);
template ifstream& operator>> (ifstream& is, CDenseArray<float>& x);
template ostream& operator<< (ostream& os, const CDenseArray<int>& x);
template istream& operator>> (istream& is, CDenseArray<int>& x);
template ifstream& operator>> (ifstream& is, CDenseArray<int>& x);
template ostream& operator<< (ostream& os, const CDenseArray<size_t>& x);
template istream& operator>> (istream& is, CDenseArray<size_t>& x);
template ifstream& operator>> (ifstream& is, CDenseArray<size_t>& x);
template ostream& operator<< (ostream& os, const CDenseArray<bool>& x);
template istream& operator>> (istream& is, CDenseArray<bool>& x);
template ifstream& operator>> (ifstream& is, CDenseArray<bool>& x);

template <class T>
CDenseVector<T>::CDenseVector():
	CDenseArray<T>::CDenseArray() {}

template <class T>
CDenseVector<T>::CDenseVector(size_t n):
	CDenseArray<T>::CDenseArray(n,1)
{}

template <class T>
CDenseVector<T>::CDenseVector(size_t nrows, size_t ncols):
    CDenseArray<T>::CDenseArray(nrows,ncols) {

    assert(nrows==1 || ncols==1);

}

template <class T>
CDenseVector<T>::CDenseVector(const CDenseVector& vector):
    CDenseArray<T>::CDenseArray(vector) {}

template <class T>
CDenseVector<T>::CDenseVector(size_t n, shared_ptr<T> data):
	CDenseArray<T>::CDenseArray(n,1,data){}

template <class T>
CDenseVector<T> CDenseVector<T>::operator=(const CDenseVector<T>& vector) {

    if(this==&vector)
        return *this;

    m_nrows = vector.m_nrows;
    m_ncols = vector.m_ncols;
    m_transpose = vector.m_transpose;
    m_data = vector.m_data;

    return *this;

}

template <class T>
CDenseVector<T> CDenseVector<T>::Clone() {

    CDenseVector<T> result(m_nrows,m_ncols);

    // copy data
    T* newdata = result.Data().get();
    T* thisdata = m_data.get();
    memcpy(newdata,thisdata,m_nrows*m_ncols*sizeof(T));

    return result;

}

template <class T>
T& CDenseVector<T>::operator()(size_t i) {

	if(m_transpose)
		return CDenseArray<T>::operator()(0,i);
	else
		return CDenseArray<T>::operator()(i,0);

}

template <class T>
T CDenseVector<T>::Get(size_t i) const {

	if(m_transpose)
		return CDenseArray<T>::Get(0,i);
	else
		return CDenseArray<T>::Get(i,0);

}

template <class T>
CDenseVector<T> CDenseVector<T>::operator+(const T& scalar) const {

    CDenseVector<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();
    T* pdatares = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = pdata[i] + scalar;

	return result;

}

template <class T>
CDenseVector<T> CDenseVector<T>::operator+(const CDenseVector<T>& vector) const {

	assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    CDenseVector<T> result(m_nrows,m_ncols);

    T* px = m_data.get();
    T* py = vector.m_data.get();
    T* pz = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pz[i] = px[i] + py[i];

	return result;

}

template <class T>
void CDenseVector<T>::Add(const CDenseVector<T>& vector) {

    assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    T* px = m_data.get();
    T* py = vector.m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        px[i] += py[i];

}

template <class T>
CDenseVector<T> CDenseVector<T>::operator-(const CDenseVector<T>& vector) const {

	assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    CDenseVector<T> result(m_nrows,m_ncols);

    T* px = m_data.get();
    T* py = vector.m_data.get();
    T* pz = result.m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pz[i] = px[i] - py[i];

	return result;

}

template <class T>
CDenseVector<T> CDenseVector<T>::operator*(const T& scalar) const {

    CDenseVector<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();
    T* pdatares = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = scalar*pdata[i];

	return result;

}

template <class T>
CDenseVector<T> CDenseVector<T>::CrossProduct(const CDenseVector<T>& x, const CDenseVector<T>& y) {

	assert((x.NRows()==3 && y.NRows()==3) || (x.NCols()==3 && y.NCols()==3));

	CDenseVector<T> result(3);

	result(0) = x.Get(1)*y.Get(2) - x.Get(2)*y.Get(1);
	result(1) = x.Get(2)*y.Get(0) - x.Get(0)*y.Get(2);
	result(2) = x.Get(0)*y.Get(1) - x.Get(1)*y.Get(0);

	return result;

}

template <class T>
void CDenseVector<T>::Sort() {

	#ifdef _OPENMP
    __gnu_parallel::sort(m_data.get(),m_data.get()+CDenseArray<T>::NElems());
	#else
    sort(m_data.get(),m_data.get()+CDenseArray<T>::NElems());
	#endif

}

template class CDenseVector<double>;
template class CDenseVector<float>;
template class CDenseVector<int>;
template class CDenseVector<size_t>;
template class CDenseVector<bool>;


template <class T>
CDenseSymmetricArray<T>::CDenseSymmetricArray():
	m_nrows(0),
	m_data(0) {

}

template <class T>
CDenseSymmetricArray<T>::CDenseSymmetricArray(size_t nrows):
	m_nrows(nrows) {

	m_data = new T[(nrows*(nrows+1))/2];

	memset(m_data,0,((nrows*(nrows+1))/2)*sizeof(T));

}

template <class T>
CDenseSymmetricArray<T>::CDenseSymmetricArray(const CDenseSymmetricArray& array):
	m_nrows(array.m_nrows) {

	m_data = new T[(array.m_nrows*(array.m_nrows+1))/2];

	memcpy(m_data,array.m_data,((array.m_nrows*(array.m_nrows+1))/2)*sizeof(T));

}

template <class T>
CDenseSymmetricArray<T>::~CDenseSymmetricArray() {

	delete [] m_data;

}

template <class T>
void CDenseSymmetricArray<T>::Print() const {

	for(size_t i=0; i<m_nrows; i++) {

		for(size_t j=0; j<m_nrows; j++) {

			printf("%.2f\t",(float)Get(i,j));

		}

		printf("\n");

	}

}

template <class T>
T CDenseSymmetricArray<T>::Norm2() const {

	T sum = 0;

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		sum += m_data[i]*m_data[i];

	return sqrt(sum);

}



template <class T>
T& CDenseSymmetricArray<T>::operator()(size_t i, size_t j) {

    assert(i<m_nrows && j<m_nrows);

	if(i>j)
		return operator()(j,i);

	return m_data[(j*(j+1))/2 + i];

}

template <class T>
T CDenseSymmetricArray<T>::Get(size_t i, size_t j) const {

    assert(i<m_nrows && j<m_nrows);

	if(i>j)
		return Get(j,i);

	return m_data[(j*(j+1))/2 + i];

}

template <class T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator=(const CDenseSymmetricArray<T>& array) {

	if(this==&array)
		return *this;

	if(m_nrows!=array.m_nrows) {

		delete [] m_data;

		m_nrows = array.m_nrows;
		m_data = new T[(array.m_nrows*(array.m_nrows+1))/2];

	}

	memcpy(m_data,array.m_data,((array.m_nrows*(array.m_nrows+1))/2)*sizeof(T));

	return *this;

}

template <class T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator+(const CDenseSymmetricArray& array) const {

	assert(m_nrows==array.m_nrows);

	CDenseSymmetricArray<T> result = CDenseSymmetricArray<T>(array.m_nrows);

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		result.m_data[i] = m_data[i] + array.m_data[i];

	return result;

}

template <class T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator-(const CDenseSymmetricArray& array) const {

	assert(m_nrows==array.m_nrows);
	CDenseSymmetricArray<T> result = CDenseSymmetricArray<T>(array.m_nrows);

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		result.m_data[i] = m_data[i] - array.m_data[i];

	return result;

}

template <class T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator*(const T& scalar) const {

	CDenseSymmetricArray<T> result = CDenseSymmetricArray<T>(m_nrows);

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		result.m_data[i] = scalar* m_data[i];

	return result;

}

template <class T>
CDenseArray<T> CDenseSymmetricArray<T>::operator*(const CDenseArray<T>& array) const {

	assert(m_nrows==array.NRows());

	CDenseArray<T> result = CDenseArray<T>(m_nrows,array.NCols());

	for(size_t i=0; i<m_nrows; i++) {

		for(size_t j=0; j<array.NCols(); j++) {

			T sum = 0;

			for(size_t k=0; k<array.NRows(); k++)
					sum += Get(i,k)*(array.Get(k,j));

			result(i,j) = sum;

		}

	}

	return result;

}

template <class T>
void CDenseSymmetricArray<T>::Scale(T scalar) {

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		m_data[i] = scalar*m_data[i];

}

template class CDenseSymmetricArray<float>;
template class CDenseSymmetricArray<double>;
template class CDenseSymmetricArray<int>;

}
