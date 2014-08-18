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

#ifdef HAVE_TBB
#include <parallel/algorithm>
#endif

#ifdef HAVE_EXR
#include <half.h>
#endif

using namespace std;

namespace R4R {

template <typename T>
CDenseArray<T>::CDenseArray():
	m_nrows(0),
	m_ncols(0),
    m_transpose(false),
    m_data() {

}

template <typename T>
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

template <typename T>
CDenseArray<T>::CDenseArray(size_t nrows, size_t ncols, const CDenseVector<T>& x):
    m_nrows(nrows),
    m_ncols(ncols),
    m_transpose(false),
    m_data(x.Data()){

    assert(nrows*ncols==x.NElems());

}

template <typename T>
void CDenseArray<T>::Resize(size_t nrows, size_t ncols) {

    if(nrows==m_nrows && ncols==m_ncols)
        return;

    CDenseArray<T> result(nrows,ncols);

    // copy data
    for(size_t i=0; i<min(nrows,m_nrows); i++) {

        for(size_t j=0; j<min(ncols,m_ncols); j++)
            result(i,j) = this->Get(i,j);

    }

    *this = result;

}

template <typename T>
void CDenseArray<T>::Concatenate(const CDenseArray& array, bool direction) {

    // put on top of each other
    if(!direction) {

        assert(m_ncols==array.NCols());

        size_t oldnrows = m_nrows;

        // now resize, this changes m_nrows
        this->Resize(oldnrows+array.NRows(),m_ncols);

        // copy
        for(size_t i=0; i<array.NRows(); i++) {

            for(size_t j=0; j<array.NCols(); j++)
                this->operator ()(oldnrows+i,j) = array.Get(i,j);

        }

    }
    else {

        assert(m_nrows==array.NRows());

        size_t oldncols = m_ncols;

        // now resize, this changes m_nrows
        this->Resize(m_nrows,oldncols+array.NCols());

        // copy
        for(size_t i=0; i<array.NRows(); i++) {

            for(size_t j=0; j<array.NCols(); j++)
                this->operator ()(i,oldncols+j) = array.Get(i,j);

        }


    }

}

template <typename T>
CDenseArray<T>::CDenseArray(const CDenseArray& array):
	m_nrows(array.m_nrows),
	m_ncols(array.m_ncols),
    m_transpose(array.m_transpose),
    m_data(array.m_data) {

}

template <typename T>
CDenseArray<T> CDenseArray<T>::operator=(const CDenseArray<T>& array) {

    if(this==&array)
        return *this;

    m_nrows = array.m_nrows;
    m_ncols = array.m_ncols;
    m_transpose = array.m_transpose;
    m_data = array.m_data;

    return *this;

}

template <typename T>
CDenseArray<T>::CDenseArray(size_t nrows, size_t ncols, shared_ptr<T> data):
	m_nrows(nrows),
	m_ncols(ncols),
    m_transpose(false),
    m_data(data) {}

template <typename T>
void CDenseArray<T>::Set(shared_ptr<T> data) {

    m_data = data;

}

template <typename T>
void CDenseArray<T>::Set(const CBivariateFunction& f) {

    for(size_t i=0; i<m_nrows; i++) {

        for(size_t j=0; j<m_ncols; j++)
            this->Set(i,j,f(j,i));

    }

}

template <typename T>
CDenseArray<T> CDenseArray<T>::Clone() const {

    CDenseArray<T> result(m_nrows,m_ncols);

    // copy data
    T* newdata = result.Data().get();
    const T* thisdata = m_data.get();
    memcpy(newdata,thisdata,m_nrows*m_ncols*sizeof(T));

    return result;

}

template <typename T>
CDenseArray<T>::~CDenseArray() {

    // nothing to do with smart pointers

}

template <typename T>
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

template <typename T>
void CDenseArray<T>::Ones() {

    T* pdata = m_data.get();

    fill_n(pdata,NElems(),1);

}

template <typename T>
void CDenseArray<T>::Zeros() {

    T* pdata = m_data.get();

    fill_n(pdata,NElems(),0);

}

template <class U>
ostream& operator<< (ostream& os, const CDenseArray<U>& x) {

	os.setf(ios::scientific,ios::floatfield);
	os.precision(5);

	for(size_t i=0; i<x.m_nrows; i++) {

		for(size_t j=0; j<x.m_ncols-1; j++) {

            os << x.Get(i,j) << " ";

		}

		if(i<x.m_nrows-1)
            os << x.Get(i,x.m_ncols-1) << endl;
		else
            os << x.Get(i,x.m_ncols-1);

	}

	return os;

}

template <>
ostream& operator<< (ostream& os, const CDenseArray<unsigned char>& x) {

    for(size_t i=0; i<x.m_nrows; i++) {

        for(size_t j=0; j<x.m_ncols-1; j++) {

            os << (unsigned int)x.Get(i,j) << " ";

        }

        if(i<x.m_nrows-1)
            os << (unsigned int)x.Get(i,x.m_ncols-1) << endl;
        else
            os << (unsigned int)x.Get(i,x.m_ncols-1);

    }

    return os;

}

template ostream& operator<< (ostream& os, const CDenseArray<double>& x);
template ostream& operator<< (ostream& os, const CDenseArray<float>& x);
template ostream& operator<< (ostream& os, const CDenseArray<int>& x);
template ostream& operator<< (ostream& os, const CDenseArray<size_t>& x);
template ostream& operator<< (ostream& os, const CDenseArray<bool>& x);
template ostream& operator<< (ostream& os, const CDenseArray<rgb>& x);
template ostream& operator<< (ostream& os, const CDenseArray<vec3>& x);

template <class U>
ofstream& operator<< (ofstream& os, const CDenseArray<U>& x) {

    os << x.NRows() << " " << x.NCols() << " " << int(GetEType<U>()) << endl;

    os.write((char*)(x.m_data.get()),sizeof(U)*x.NElems());

    return os;

}

template ofstream& operator<< (ofstream& os, const CDenseArray<double>& x);
template ofstream& operator<< (ofstream& os, const CDenseArray<float>& x);
template ofstream& operator<< (ofstream& os, const CDenseArray<int>& x);
template ofstream& operator<< (ofstream& os, const CDenseArray<size_t>& x);
template ofstream& operator<< (ofstream& os, const CDenseArray<bool>& x);
template ofstream& operator<< (ofstream& os, const CDenseArray<rgb>& x);
template ofstream& operator<< (ofstream& os, const CDenseArray<vec3>& x);

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

    int temp;
    in >> temp;
    ETYPE type = (ETYPE)temp;
    in.get();

    U* pdata = x.m_data.get();

    size_t nelems = nrows*ncols;

    switch(type) {

    case ETYPE::B1U:
    {

        bool* buffer = new bool[nelems];
        in.read((char*)buffer,nelems*sizeof(bool));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::C1S:
    {

        char* buffer = new char[nelems];
        in.read((char*)buffer,nelems*sizeof(char));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::C1U:
    {

        unsigned char* buffer = new unsigned char[nelems];
        in.read((char*)buffer,nelems*sizeof(unsigned char));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::S2S:
    {

        short* buffer = new short[nelems];
        in.read((char*)buffer,nelems*sizeof(short));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }
    case ETYPE::S2U:
    {

        unsigned short* buffer = new unsigned short[nelems];
        in.read((char*)buffer,nelems*sizeof(unsigned short));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::I4S:
    {

        int* buffer = new int[nelems];
        in.read((char*)buffer,nelems*sizeof(int));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::I4U:
    {

        unsigned int* buffer = new unsigned int[nelems];
        in.read((char*)buffer,nelems*sizeof(unsigned int));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::F4S:
    {

        float* buffer = new float[nelems];
        in.read((char*)buffer,nelems*sizeof(float));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::L8S:
    {

        long int* buffer = new long int[nelems];
        in.read((char*)buffer,nelems*sizeof(long int));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::L8U:
    {

        size_t* buffer = new size_t[nelems];
        in.read((char*)buffer,nelems*sizeof(size_t));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    case ETYPE::D8S:
    {

        double* buffer = new double[nelems];
        in.read((char*)buffer,nelems*sizeof(double));

        // copy data and cast
        for(size_t i=0; i<nelems; i++)
            pdata[i] = (U)buffer[i];

        delete [] buffer;

        break;

    }

    default:
        return in;


    }

    return in;

}

template ifstream& operator >> (ifstream& is, CDenseArray<double>& x);
template ifstream& operator >> (ifstream& is, CDenseArray<float>& x);
template ifstream& operator >> (ifstream& is, CDenseArray<int>& x);
template ifstream& operator >> (ifstream& is, CDenseArray<size_t>& x);
template ifstream& operator >> (ifstream& is, CDenseArray<bool>& x);

template <>
ifstream& operator >> (ifstream& in, CDenseArray<vec3>& x) {

    // read information
    size_t nrows, ncols;
    in >> nrows;
    in >> ncols;

    // resize storage if necessary
    if(nrows!=x.NRows() || ncols!=x.NCols()) {

        x.m_data.reset();

#ifndef __SSE4_1__
        x.m_data.reset(new vec3[nrows*ncols]);
#else
        x.m_data.reset((vec3*)_mm_malloc(nrows*ncols*sizeof(vec3),16));
#endif
        x.m_nrows = nrows;
        x.m_ncols = ncols;
        x.m_transpose = false;

    }

    int temp;
    in >> temp;
    ETYPE type = (ETYPE)temp;
    in.get();

    if(type!=ETYPE::D8S3) {

        cerr << "Could not read file..." << endl;

    }

    in.read((char*)x.m_data.get(),nrows*ncols*sizeof(vec3));

    return in;

}

template ifstream& operator >> (ifstream& is, CDenseArray<vec3>& x);

template <typename T>
bool CDenseArray<T>::WriteToFile(const char* filename) {

    ofstream out(filename);

    if(!out.is_open()) {

        cerr << "ERROR: Could not open " << filename << "..." << endl;
        return 1;

    }

    out << *this;

    out.close();

    return 0;

}

template <typename T>
bool CDenseArray<T>::ReadFromFile(const char* filename) {

    ifstream in(filename);

    if(!in.is_open()) {

        cerr << "ERROR: File " << filename << " not found..." << endl;
        return 1;

    }

    in >> *this;

    in.close();

    return 0;

}

template <typename T>
bool CDenseArray<T>::Normalize() {

    double norm = Norm2();

	if(norm>0) {

        T s = T(1/norm);

		Scale(s);

		return 0;

	}

	return 1;

}


template <typename T>
T CDenseArray<T>::Trace() const {

	size_t m = min(m_nrows,m_ncols);

	T sum = 0;

	for(size_t i=0; i<m; i++)
		sum += Get(i,i);

	return sum;

}

template <typename T>
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

template <typename T>
bool CDenseArray<T>::Invert() {

    // this is not implemented yet for fields other than the reals
    assert(false);

    return false;

}

template <>
bool CDenseArray<double>::Invert() {

    assert((m_ncols==2 && m_nrows==2) || (m_ncols==3 && m_nrows==3));

    double* pdata = m_data.get();
    double* temp = new double[m_nrows*m_ncols];
    memcpy(temp,pdata,m_nrows*m_ncols*sizeof(double));

    if(m_nrows==2) {

        double det = temp[0]*temp[3]-temp[1]*temp[2];

        if(det==0) {

            delete [] temp;
            return 1;

        }

        double deti = 1/det;

        // first linear part
        pdata[0] = deti*temp[3];
        pdata[1] = -deti*temp[1];
        pdata[2] = -deti*temp[2];
        pdata[3] = deti*temp[0];

    }
    else if(m_nrows==3) {

        // inverse of det
        double det = temp[0]*(temp[4]*temp[8] - temp[5]*temp[7])
                   - temp[3]*(temp[8]*temp[1] - temp[2]*temp[7])
                   + temp[6]*(temp[1]*temp[5] - temp[2]*temp[4]);

        if(det==0) {

            delete [] temp;
            return 1;

        }

        double deti = 1.0/det;

        // first linear part
        pdata[0] = deti*(temp[4]*temp[8] - temp[5]*temp[7]);
        pdata[1] = -deti*(temp[1]*temp[8] - temp[2]*temp[7]);
        pdata[2] = deti*(temp[1]*temp[5] - temp[2]*temp[4]);
        pdata[3] = -deti*(temp[3]*temp[8] - temp[5]*temp[6]);
        pdata[4] = deti*(temp[0]*temp[8] - temp[2]*temp[6]);
        pdata[5] = -deti*(temp[0]*temp[5] - temp[2]*temp[3]);
        pdata[6] = deti*(temp[3]*temp[7] - temp[4]*temp[6]);
        pdata[7] = -deti*(temp[0]*temp[7] - temp[1]*temp[6]);
        pdata[8] = deti*(temp[0]*temp[4] - temp[1]*temp[3]);

    }
    else {

        delete [] temp;
        return 1;

    }

    delete [] temp;

    return 0;

}

template <typename T>
void CDenseArray<T>::Transpose() {

    std::swap(m_nrows,m_ncols);
	m_transpose = !m_transpose;

}

template <typename T>
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

template <typename T>
double CDenseArray<T>::Norm1() const {

    CHellingerKernel<T> kernel(m_nrows*m_ncols);

    T* pdata = m_data.get();

    return kernel.ComputeKernelNorm(pdata);

}

template <typename T>
double CDenseArray<T>::Norm(double p) const {

    assert(p>0);

	double sum = 0;

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        sum += pow(fabs(pdata[i]),p);

    return pow(sum,1/p);

}

template <>
double CDenseArray<vec3>::Norm(double p) const {

    assert(p>0);

    double sum = 0;

    vec3* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++)
            sum += pow(fabs(pdata[i].Get(j)),p);

    }

    return pow(sum,1/p);

}

template <>
double CDenseArray<vec3f>::Norm(double p) const {

    assert(p>0);

    double sum = 0;

    vec3f* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++)
            sum += pow(fabs(pdata[i].Get(j)),p);

    }

    return pow(sum,1/p);

}

template <>
double CDenseArray<rgb>::Norm(double p) const {

    assert(p>0);

    double sum = 0;

    rgb* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++)
            sum += pow(fabs(pdata[i].Get(j)),p);

    }

    return pow(sum,1/p);

}

template <typename T>
T CDenseArray<T>::Sum() const {

	T sum = 0;

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        sum += pdata[i];

	return sum;

}


template <typename T>
void CDenseArray<T>::Abs() {

    T* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = fabs(pdata[i]);

}

template <>
void CDenseArray<vec3>::Abs() {

    vec3* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = pdata[i].Abs();

}


template <>
void CDenseArray<vec3f>::Abs() {

    vec3f* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = pdata[i].Abs();

}
template <>
void CDenseArray<rgb>::Abs() {

    // unsigned char are positive already, so do nothing
    return;

}

template <typename T>
T CDenseArray<T>::Get(size_t i, size_t j) const {

    assert(i<m_nrows && j<m_ncols);

    T* pdata = m_data.get();

	if(!m_transpose)
        return pdata[m_nrows*j + i];
	else
        return pdata[m_ncols*i + j];

}

/* notes:
 * - we cannot specialize this template for integral types because we
 * would have to specialize the entire class.
 * - floating point inputs are implicitly cast to doubles.
 */
template <typename T>
template<typename U>
U CDenseArray<T>::Get(const CVector<double,2> &p) const {

    double u, v;
    u = p.Get(0);
    v = p.Get(1);

    int i = (int)floor(v);
    int j = (int)floor(u);

    if(i<0 || i>NRows()-2 || j<0 || j>NCols()-2)
        return U(0);

    double vd = v - i;
    double ud = u - j;

    U A00, A01, A10, A11, A0, A1;
    A00 = U(Get(i,j));
    A01 = U(Get(i,j+1));
    A10 = U(Get(i+1,j));
    A11 = U(Get(i+1,j+1));

    A0 = A00*(1-ud) + A01*ud;
    A1 = A10*(1-ud) + A11*ud;

    return A0*(1-vd) + A1*vd;

}

template double CDenseArray<double>::Get<double>(const CVector<double,2>& p) const;
template double CDenseArray<float>::Get<double>(const CVector<double,2>& p) const;
template float CDenseArray<float>::Get<float>(const CVector<double,2>& p) const;
template double CDenseArray<unsigned char>::Get<double>(const CVector<double,2>& p) const;
template vec3 CDenseArray<rgb>::Get<vec3>(const CVector<double,2>& p) const;
template vec3 CDenseArray<vec3>::Get<vec3>(const CVector<double,2>& p) const;

template <typename T>
template<typename U>
vector<U> CDenseArray<T>::Gradient(const CVector<double,2>& p) const {

    if(p.Get(0)<1 || p.Get(1)<1 || p.Get(0)>=NCols()-2 || p.Get(1)>=NRows()-2)
        return vector<U>(2);

    vec2 dx = { 1, 0 };
    vec2 dy = { 0, 1 };

    U I0x, I0y, I1x, I1y;
    I0x = this->Get<U>(p-dx);
    I0y = this->Get<U>(p-dy);
    I1x = this->Get<U>(p+dx);
    I1y = this->Get<U>(p+dy);

    vector<U> grad;
    grad.push_back(0.5*(I1x-I0x));
    grad.push_back(0.5*(I1y-I0y));

    return grad;

}

template vector<double> CDenseArray<double>::Gradient<double>(const CVector<double,2>& p) const;
template vector<double> CDenseArray<unsigned char>::Gradient<double>(const CVector<double,2>& p) const;
template vector<vec3> CDenseArray<rgb>::Gradient<vec3>(const CVector<double,2>& p) const;

template <typename T>
template<typename U>
vector<U> CDenseArray<T>::Gradient(size_t i, size_t j) const {

    if(i<1 || j<1 || j>NCols()-2 || i>NRows()-2)
        return vector<U>(2);

    U I0x, I0y, I1x, I1y;
    I0x = U(this->Get(i,j-1));
    I0y = U(this->Get(i-1,j));
    I1x = U(this->Get(i,j+1));
    I1y = U(this->Get(i+1,j));

    vector<U> grad;
    grad.push_back(0.5*(I1x-I0x));
    grad.push_back(0.5*(I1y-I0y));

    return grad;
}

template vector<double> CDenseArray<double>::Gradient<double>(size_t i, size_t j) const;
template vector<vec3> CDenseArray<rgb>::Gradient<vec3>(size_t i, size_t j) const;

//template <typename T>
//template <typename U>
//CVector<U,2> CDenseArray<T>::ProjectToBoundary(const CVector<U,2>& x) const {

//    CVector<U,2> result = x;

//    if(x.Get(0)<0)
//        result(0) = 0;

//    if(x.Get(0)>=this->NCols())
//        result(0) = this->NCols() - 1;

//    if(x.Get(1)<0)
//        result(1) = 0;

//    if(x.Get(1)>=this->NRows())
//        result(1) = this->NRows() - 1;

//    return result;

//}

template <typename T>
CDenseVector<T> CDenseArray<T>::GetColumn(size_t j) const {

    if(m_transpose)
        return GetRow(j);

    // if sizeof(T)*nrows is a multiple of 16, we can make a shallow copy
    if(sizeof(T)*m_nrows%16==0) {

        shared_ptr<T> colptr(m_data,m_data.get()+j*m_nrows);

        // create column view
        CDenseVector<T> col(m_nrows,colptr);

        return col;

    }

    CDenseVector<T> col(m_nrows); // this will be 16-byte aligned

    for(size_t i=0; i<m_nrows; i++)
        col(i) = Get(i,j);

    return col;

}

template <typename T>
void CDenseArray<T>::SetColumn(size_t j, const CDenseVector<T>& col) {

    assert(col.NCols()==1 && col.NRows() == NRows() && j<=NCols());

	for(size_t i=0; i<col.NRows(); i++)
		this->operator ()(i,j) = col.Get(i);

}

template <typename T>
void CDenseArray<T>::SetRow(size_t i, const CDenseVector<T>& row) {

    assert(row.NRows()==1 && row.NCols() == NCols() && i<=NRows());

    for(size_t j=0; j<row.NCols(); j++)
        this->operator ()(i,j) = row.Get(j);

}

template <typename T>
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

template <typename T>
T& CDenseArray<T>::operator()(size_t i, size_t j) {

    assert(i<m_nrows && j<m_ncols);

    T* pdata = m_data.get();

	if(!m_transpose)
        return pdata[m_nrows*j + i];
	else
        return pdata[m_nrows*i + j];				// m_nrows is replaced by m_ncols in transpose method

}

template <typename T>
CDenseArray<T> CDenseArray<T>::operator+(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = pdata[i] + scalar;

	return result;

}

template <typename T>
CDenseArray<T> CDenseArray<T>::operator/(const CDenseArray<T>& array) const {

    assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

    CDenseArray<T> result(array.m_nrows,array.m_ncols);

    // this does not work low-level because of possible transpositions
    for(size_t i=0; i<m_nrows; i++) {

        for(size_t j=0; j<m_ncols; j++) {

            result(i,j) = this->Get(i,j)/array.Get(i,j);

        }

    }

    return result;

}

template <typename T>
CDenseArray<T> CDenseArray<T>::operator-(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();
    T* pdatares = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = pdata[i] - scalar;

	return result;

}


template <typename T>
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

template <typename T>
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

template <typename T>
double CDenseArray<T>::InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y) {

	assert(x.m_nrows==y.m_nrows && x.m_ncols==y.m_ncols);

    T* px = x.m_data.get();
    T* py = y.m_data.get();

    CMercerKernel<T> kernel(x.m_nrows*x.m_ncols);

    return kernel.Evaluate(px,py);

}

template <typename T>
double CDenseArray<T>::InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y, CMercerKernel<T>& kernel) {

    assert(x.m_nrows==y.m_nrows && x.m_ncols==y.m_ncols);

    T* px = x.m_data.get();
    T* py = y.m_data.get();

    return kernel.Evaluate(px,py);

}

template <typename T>
CDenseVector<T> CDenseArray<T>::ColumwiseInnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y) {

    assert(x.m_nrows==y.m_nrows && x.m_ncols==y.m_ncols);

    CDenseVector<T> result(x.m_ncols);

    for(size_t i=0; i<x.m_nrows; i++) {

        for(size_t j=0; j<x.m_ncols; j++)
            result(j) += x.Get(i,j)*y.Get(i,j);

    }

    return result;

}


template <typename T>
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

template <typename T>
CDenseArray<T> CDenseArray<T>::Transpose(const CDenseArray<T>& x) {

	CDenseArray<T> result(x);

	result.Transpose();

	return result;

}


template <typename T>
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

template <typename T>
CDenseArray<T> CDenseArray<T>::operator*(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdatares = result.m_data.get();
    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = scalar*pdata[i];

	return result;

}

template <typename T>
CDenseArray<T> CDenseArray<T>::operator/(const T& scalar) const {

    CDenseArray<T> result(m_nrows,m_ncols);

    T* pdatares = result.m_data.get();
    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = pdata[i]/scalar;

	return result;

}

template <typename T>
CDenseArray<T> CDenseArray<T>::operator*(const CDenseArray<T>& array) const {

	assert(m_ncols==array.m_nrows);

    CDenseArray<T> result(m_nrows,array.m_ncols);

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

template <typename T>
CDenseVector<T> CDenseArray<T>::operator*(const CDenseVector<T>& vector) const {

	assert(m_ncols==vector.m_nrows);

    CDenseVector<T> result(m_nrows);

	for(size_t i=0; i<m_nrows; i++) {

		T sum = 0;

		for(size_t k=0; k<m_ncols; k++)
			sum += Get(i,k)*(vector.Get(k));

		result(i) = sum;

	}

	return result;

}

/*template<typename T>
template<class Array> Array CDenseArray<T>::operator*(const Array& array) const {


    assert(m_ncols==array.m_nrows);

    Array result(m_nrows,array.m_ncols);

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

template CDenseArray<double> CDenseArray<double>::operator *(const CDenseArray<double>& x) const;
template CDenseVector<double> CDenseArray<double>::operator *(const CDenseVector<double>& x) const;
template CDenseArray<float> CDenseArray<float>::operator *(const CDenseArray<float>& x) const;
template CDenseVector<float> CDenseArray<float>::operator *(const CDenseVector<float>& x) const;
template CDenseVector<int> CDenseArray<int>::operator *(const CDenseVector<int>& x) const;
template CDenseArray<int> CDenseArray<int>::operator *(const CDenseArray<int>& x) const;*/


template<typename T>
template <u_int n> CVector<T,n> CDenseArray<T>::operator*(const CVector<T,n>& vector) const {

    assert(m_ncols==n && m_nrows==n);

    CVector<T,n> result;

    for(size_t i=0; i<n; i++) {

        T sum = 0;

        for(size_t k=0; k<m_ncols; k++)
            sum += Get(i,k)*(vector.Get(k));

        result(i) = sum;

    }

    return result;

}

template CVector<double,2> CDenseArray<double>::operator*(const CVector<double,2>& vector) const;
template CVector<double,3> CDenseArray<double>::operator*(const CVector<double,3>& vector) const;
template CVector<float,2> CDenseArray<float>::operator*(const CVector<float,2>& vector) const;
template CVector<float,3> CDenseArray<float>::operator*(const CVector<float,3>& vector) const;

template <typename T>
void CDenseArray<T>::Scale(T scalar) {

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] *= scalar;

}

template <typename T>
CDenseArray<T> CDenseArray<T>::ScaleColumns(const CDenseVector<T>& s) {

    CDenseArray<T> result(m_nrows,m_ncols);

    for(size_t i=0; i<m_nrows; i++) {

        for(size_t j=0; j<m_ncols; j++)
            result(i,j) = this->Get(i,j)*s.Get(j);

    }

    return result;

}

template <typename T>
void CDenseArray<T>::Add(const T& scalar) {

    T* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] += scalar;

}

template <typename T>
void CDenseArray<T>::Subtract(const T& scalar) {

    T* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = pdata[i] - scalar;

}

template <typename T>
T CDenseArray<T>::Median() const {

    CDenseArray<T> temp = this->Clone();

#ifdef HAVE_TBB
    __gnu_parallel::sort(temp.m_data.get(),temp.m_data.get()+temp.NElems());
#else
    sort(temp.m_data.get(),temp.m_data.get()+temp.NElems());
#endif

	if(temp.NElems()%2==1)
        return temp.m_data.get()[(size_t)((temp.NElems()+1)/2)-1];
	else
        return 0.5*(temp.m_data.get()[temp.NElems()/2-1]+temp.m_data.get()[temp.NElems()/2]);

}

template <typename T>
T CDenseArray<T>::Variance() const {

	T mean = this->Mean();

    CDenseArray<T> temp = this->Clone();

    T* pdata = temp.m_data.get();
    const T* pthisdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = (pthisdata[i] - mean)*(pthisdata[i] - mean);

	return temp.Mean();

}

template <typename T>
T CDenseArray<T>::MAD() const {

    T median = Median();

    CDenseArray<T> temp = this->Clone();

    T* pdata = temp.m_data.get();
    const T* pthisdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdata[i] = fabs(pthisdata[i]-median);

    return temp.Median();

}

template <>
rgb CDenseArray<rgb>::MAD() const {

    rgb median = Median();

    CDenseArray<rgb> temp = this->Clone();

    rgb* pdata = temp.m_data.get();
    const rgb* pthisdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++)
            pdata[i](j) = (unsigned char)fabs(pthisdata[i].Get(j)-median.Get(j));

    }

    return temp.Median();

}

template <>
vec3 CDenseArray<vec3>::MAD() const {

    //! FIXME: This is not per channel but w.r.t. to order-relation of vec3.
    vec3 median = Median();

    CDenseArray<vec3> temp = this->Clone();

    vec3* pdata = temp.m_data.get();
    const vec3* pthisdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++)
            pdata[i](j) = fabs(pthisdata[i].Get(j)-median.Get(j));

    }

    return temp.Median();

}

template <>
vec3f CDenseArray<vec3f>::MAD() const {

    //! FIXME: This is not per channel but w.r.t. to order-relation of vec3.
    vec3f median = Median();

    CDenseArray<vec3f> temp = this->Clone();

    vec3f* pdata = temp.m_data.get();
    const vec3f* pthisdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++)
            pdata[i](j) = fabs(pthisdata[i].Get(j)-median.Get(j));

    }

    return temp.Median();

}

template <typename T>
T CDenseArray<T>::Min() const {

	T min = numeric_limits<T>::max();

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++) {

        if(pdata[i]<=min)
            min = pdata[i];

	}

	return min;

}

template <>
rgb CDenseArray<rgb>::Min() const {

    rgb min = { 0, 0, 0 };
    rgb* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++) {

            // find the minimum of every channel
            if(pdata[i].Get(j)<=min.Get(j))
                min(j) = pdata[i].Get(j);

        }

    }

    return min;

}

template <typename T>
T CDenseArray<T>::Max() const {

    T max = numeric_limits<T>::max()*(-1);

    T* pdata = m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++) {

        if(max<=pdata[i])
            max = pdata[i];

	}

	return max;

}

template <>
rgb CDenseArray<rgb>::Max() const {

    rgb max = { 255, 255, 255 };
    rgb* pdata = m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++) {

        for(size_t j=0; j<3; j++) {

            // find the minimum of every channel
            if(pdata[i].Get(j)>=max.Get(j))
                max(j) = pdata[i].Get(j);

        }

    }

    return max;

}

template class CDenseArray<float>;
template class CDenseArray<double>;
template class CDenseArray<int>;
template class CDenseArray<size_t>;
template class CDenseArray<bool>;
template class CDenseArray<rgb>;
template class CDenseArray<vec3>;
template class CDenseArray<vec3f>;
template class CDenseArray<unsigned char>;

#ifdef HAVE_EXR
template class CDenseArray<half>;
//template class CDenseArray<vec3h>;
#endif

template <typename T>
CDenseVector<T>::CDenseVector():
	CDenseArray<T>::CDenseArray() {}

template <typename T>
CDenseVector<T>::CDenseVector(size_t n):
	CDenseArray<T>::CDenseArray(n,1)
{}

template <typename T>
CDenseVector<T>::CDenseVector(size_t nrows, size_t ncols):
    CDenseArray<T>::CDenseArray(nrows,ncols) {

    assert(nrows==1 || ncols==1);

    if(nrows==1)
        m_transpose = true;

}

template <typename T>
CDenseVector<T>::CDenseVector(const CDenseVector& vector):
    CDenseArray<T>::CDenseArray(vector) {}

template <typename T>
CDenseVector<T>::CDenseVector(size_t n, shared_ptr<T> data):
	CDenseArray<T>::CDenseArray(n,1,data){}

template <typename T>
CDenseVector<T> CDenseVector<T>::operator=(const CDenseVector<T>& vector) {

    if(this==&vector)
        return *this;

    m_nrows = vector.m_nrows;
    m_ncols = vector.m_ncols;
    m_transpose = vector.m_transpose;
    m_data = vector.m_data;

    return *this;

}

template <typename T>
template<u_int n>
CDenseVector<T>::CDenseVector(CVector<T,n>& x):
    CDenseArray<T>::CDenseArray(n,1) {

    memcpy(m_data.get(),x.Data(),n*sizeof(T));

}

template <typename T>
CDenseVector<T>::CDenseVector(const CDenseArray<T>& x):
    CDenseArray<T>::CDenseArray(x.NElems(),1,x.Data()){}

template CDenseVector<double>::CDenseVector<3>(CVector<double,3>& x);
template CDenseVector<double>::CDenseVector<2>(CVector<double,2>& x);
template CDenseVector<float>::CDenseVector<3>(CVector<float,3>& x);
template CDenseVector<float>::CDenseVector<2>(CVector<float,2>& x);

template <typename T>
CDenseVector<T> CDenseVector<T>::Clone() const {

    CDenseVector<T> result(m_nrows,m_ncols);

    // copy data
    T* newdata = result.Data().get();
    const T* thisdata = m_data.get();
    memcpy(newdata,thisdata,m_nrows*m_ncols*sizeof(T));

    return result;

}

template <typename T>
T& CDenseVector<T>::operator()(size_t i) {

	if(m_transpose)
		return CDenseArray<T>::operator()(0,i);
	else
		return CDenseArray<T>::operator()(i,0);

}

template <typename T>
T CDenseVector<T>::Get(size_t i) const {

	if(m_transpose)
		return CDenseArray<T>::Get(0,i);
	else
		return CDenseArray<T>::Get(i,0);

}

template <typename T>
CDenseVector<T> CDenseVector<T>::operator+(const T& scalar) const {

    CDenseVector<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();
    T* pdatares = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = pdata[i] + scalar;

	return result;

}

template <typename T>
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

template <typename T>
void CDenseVector<T>::Add(const CDenseVector<T>& vector) {

    assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    T* px = m_data.get();
    T* py = vector.m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        px[i] += py[i];

}

template <>
CDenseVector<float> CDenseVector<float>::operator+(const CDenseVector<float>& vector) const {

    assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    CDenseVector<float> result(m_nrows,m_ncols);

    float* px = m_data.get();
    float* py = vector.m_data.get();
    float* pz = result.m_data.get();

//#ifndef __SSE4_1__
    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pz[i] = px[i] + py[i];
/*#else
    size_t n = m_nrows*m_ncols;

    size_t offset = n - n%4;    // #m_n - #m_n%4

    __m128* pxs = (__m128*)px;
    __m128* pys = (__m128*)py;
    __m128* pzs = (__m128*)pz;

    for(size_t i=0; i<offset/4; i++)
        pzs[i] = _mm_hadd_ps(pxs[i],pys[i]);


    // add offset
    for(size_t i=offset; i<n; i++)
        pz[i] = px[i] + py[i];
#endif*/

    return result;

}

template <>
CDenseVector<float> CDenseVector<float>::operator-(const CDenseVector<float>& vector) const {

    assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    CDenseVector<float> result(m_nrows,m_ncols);

    float* px = m_data.get();
    float* py = vector.m_data.get();
    float* pz = result.m_data.get();

//#ifndef __SSE4_1__
    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pz[i] = px[i] - py[i];
/*#else
    size_t n = m_nrows*m_ncols;

    size_t offset = n - n%4;    // #m_n - #m_n%4

    __m128* pxs = (__m128*)px;
    __m128* pys = (__m128*)py;
    __m128* pzs = (__m128*)pz;

    for(size_t i=0; i<offset/4; i++)
        pzs[i] = _mm_hsub_ps(pxs[i],pys[i]);


    // add offset
    for(size_t i=offset; i<n; i++)
        pz[i] = px[i] - py[i];
#endif*/

    return result;

}

template <typename T>
CDenseVector<T> CDenseVector<T>::operator/(const CDenseVector<T>& vector) const {

    assert(m_nrows==vector.m_nrows && m_ncols==vector.m_ncols);

    CDenseVector<T> result(m_nrows,m_ncols);

    const T* px = m_data.get();
    const T* py = vector.m_data.get();
    T* pz = result.m_data.get();

    for(size_t i=0; i<m_nrows*m_ncols; i++)
        pz[i] = px[i]*(1/py[i]);

    return result;

}


template <typename T>
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

template <typename T>
CDenseVector<T> CDenseVector<T>::operator*(const T& scalar) const {

    CDenseVector<T> result(m_nrows,m_ncols);

    T* pdata = m_data.get();
    T* pdatares = result.m_data.get();

	for(size_t i=0; i<m_nrows*m_ncols; i++)
        pdatares[i] = scalar*pdata[i];

	return result;

}

template <typename T>
CDenseVector<T> CDenseVector<T>::CrossProduct(const CDenseVector<T>& x, const CDenseVector<T>& y) {

	assert((x.NRows()==3 && y.NRows()==3) || (x.NCols()==3 && y.NCols()==3));

	CDenseVector<T> result(3);

	result(0) = x.Get(1)*y.Get(2) - x.Get(2)*y.Get(1);
	result(1) = x.Get(2)*y.Get(0) - x.Get(0)*y.Get(2);
	result(2) = x.Get(0)*y.Get(1) - x.Get(1)*y.Get(0);

	return result;

}

template <typename T>
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
template class CDenseVector<u_char>;


template <typename T>
CDenseSymmetricArray<T>::CDenseSymmetricArray():
	m_nrows(0),
	m_data(0) {

}

template <typename T>
CDenseSymmetricArray<T>::CDenseSymmetricArray(size_t nrows):
	m_nrows(nrows) {

	m_data = new T[(nrows*(nrows+1))/2];

	memset(m_data,0,((nrows*(nrows+1))/2)*sizeof(T));

}

template <typename T>
CDenseSymmetricArray<T>::CDenseSymmetricArray(const CDenseSymmetricArray& array):
	m_nrows(array.m_nrows) {

	m_data = new T[(array.m_nrows*(array.m_nrows+1))/2];

	memcpy(m_data,array.m_data,((array.m_nrows*(array.m_nrows+1))/2)*sizeof(T));

}

template <typename T>
CDenseSymmetricArray<T>::~CDenseSymmetricArray() {

	delete [] m_data;

}

template <typename T>
void CDenseSymmetricArray<T>::Print() const {

	for(size_t i=0; i<m_nrows; i++) {

		for(size_t j=0; j<m_nrows; j++) {

			printf("%.2f\t",(float)Get(i,j));

		}

		printf("\n");

	}

}

template <typename T>
T CDenseSymmetricArray<T>::Norm2() const {

	T sum = 0;

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		sum += m_data[i]*m_data[i];

	return sqrt(sum);

}



template <typename T>
T& CDenseSymmetricArray<T>::operator()(size_t i, size_t j) {

    assert(i<m_nrows && j<m_nrows);

	if(i>j)
		return operator()(j,i);

	return m_data[(j*(j+1))/2 + i];

}

template <typename T>
T CDenseSymmetricArray<T>::Get(size_t i, size_t j) const {

    assert(i<m_nrows && j<m_nrows);

	if(i>j)
		return Get(j,i);

	return m_data[(j*(j+1))/2 + i];

}

template <typename T>
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

template <typename T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator+(const CDenseSymmetricArray& array) const {

	assert(m_nrows==array.m_nrows);

	CDenseSymmetricArray<T> result = CDenseSymmetricArray<T>(array.m_nrows);

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		result.m_data[i] = m_data[i] + array.m_data[i];

	return result;

}

template <typename T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator-(const CDenseSymmetricArray& array) const {

	assert(m_nrows==array.m_nrows);
	CDenseSymmetricArray<T> result = CDenseSymmetricArray<T>(array.m_nrows);

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		result.m_data[i] = m_data[i] - array.m_data[i];

	return result;

}

template <typename T>
CDenseSymmetricArray<T> CDenseSymmetricArray<T>::operator*(const T& scalar) const {

	CDenseSymmetricArray<T> result = CDenseSymmetricArray<T>(m_nrows);

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		result.m_data[i] = scalar* m_data[i];

	return result;

}

template <typename T>
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

template <typename T>
void CDenseSymmetricArray<T>::Scale(T scalar) {

	for(size_t i=0; i<(m_nrows*(m_nrows+1))/2; i++)
		m_data[i] = scalar*m_data[i];

}

template class CDenseSymmetricArray<float>;
template class CDenseSymmetricArray<double>;
template class CDenseSymmetricArray<int>;

}
