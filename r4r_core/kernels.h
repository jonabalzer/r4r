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

#ifndef R4RKERNELS_H
#define R4RKERNELS_H

#include <iostream>
#include <math.h>

namespace R4R {

enum KERNEL {  IDENTITY, CHISQUARED, INTERSECTION, HELLINGER };

template <class T>
class CMercerKernel {

public:

    //! Constructor.
    CMercerKernel():m_n(0), m_offset(0) {}

    //! Constructor.
    CMercerKernel(int n):m_n(n),m_offset(n-n%4){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y);

    //! Compute the norm induced by the kernel.
    double ComputeKernelNorm(T* x) { return sqrt(Evaluate(x,x)); }

    //! Computes the gradient the kernel w.r.t. \f$x\f$.
    void Gradient(T* x, T* y, T* nablax);

    //! Create a kernel by number. Make sure to de-allocate it later.
    static CMercerKernel<T>* Create(int no, int n);

    //! Tests the kernel.
    static void TestKernel(int kn, int n, size_t notests);

    //! Access to the size.
    int GetN() { return m_n; }

protected:

    int m_n;            // dimension of the target space
    int m_offset;       // #m_n - #m_n%4

};

template <class T>
class CChiSquaredKernel: public CMercerKernel<T> {

public:

    //! Constructor.
    CChiSquaredKernel():CMercerKernel<T>::CMercerKernel() {}

    //! Constructor.
    CChiSquaredKernel(int n):CMercerKernel<T>::CMercerKernel(n){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y);

private:

    using CMercerKernel<T>::m_n;
    using CMercerKernel<T>::m_offset;

};


template <class T>
class CIntersectionKernel: public CMercerKernel<T> {

public:

    //! Constructor.
    CIntersectionKernel():CMercerKernel<T>::CMercerKernel() {}

    //! Constructor.
    CIntersectionKernel(int n):CMercerKernel<T>::CMercerKernel(n){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y);

private:

    using CMercerKernel<T>::m_n;
    using CMercerKernel<T>::m_offset;

};


template <class T>
class CHellingerKernel: public CMercerKernel<T> {

public:

    //! Constructor.
    CHellingerKernel():CMercerKernel<T>::CMercerKernel() {}

    //! Constructor.
    CHellingerKernel(int n):CMercerKernel<T>::CMercerKernel(n){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y);

private:

   using CMercerKernel<T>::m_n;

};

template <class T>
class CRBFKernel: public CMercerKernel<T> {

public:

    //! Constructor.
    CRBFKernel():CMercerKernel<T>::CMercerKernel(),m_sigma(1){}

    //! Constructor.
    CRBFKernel(int n):CMercerKernel<T>::CMercerKernel(n),m_sigma(1){}

    //! Constructor.
    CRBFKernel(int n, double sigma):CMercerKernel<T>::CMercerKernel(n),m_sigma(sigma){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y);

private:

   using CMercerKernel<T>::m_n;
   double m_sigma;

};

template <class T>
class CPolynomialKernel: public CMercerKernel<T> {

public:

    //! Constructor.
    CPolynomialKernel():CMercerKernel<T>::CMercerKernel(),m_c(0),m_d(1){}

    //! Constructor.
    CPolynomialKernel(int n):CMercerKernel<T>::CMercerKernel(n),m_c(0),m_d(1){}

    //! Constructor.
    CPolynomialKernel(int n, double c, double d):CMercerKernel<T>::CMercerKernel(n),m_c(c), m_d(d){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y);

private:

   using CMercerKernel<T>::m_n;
   double m_c;
   double m_d;

};


template<class T>
class CHammingKernel: public CMercerKernel<T> {

public:

    //! Constructor.
    CHammingKernel():CMercerKernel<T>::CMercerKernel() {}

    //! Constructor.
    CHammingKernel(int n):CMercerKernel<T>::CMercerKernel(n){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y) {

        size_t result = 0;

        for(int i=0; i<m_n; i++) {

            if(x[i]!=y[i])
                (i+1)*result++;

        }

        return (double)(m_n-result);

    }

private:

    using CMercerKernel<T>::m_n;

};

//! Interprets the input as bitstreams and computes their distance.
template<class T>
class CBitwiseHammingKernel: public CMercerKernel<T> {

public:

    static unsigned int m_lut[256];

    //! Constructor.
    CBitwiseHammingKernel():CMercerKernel<T>::CMercerKernel() {}

    //! Constructor.
    CBitwiseHammingKernel(int n):CMercerKernel<T>::CMercerKernel(n){}

    //! Evaluate kernel.
    virtual double Evaluate(T* x, T* y) {

        size_t result = 0;

        unsigned char* px = (unsigned char*)x;
        unsigned char* py = (unsigned char*)y;

        for(int i=0; i<m_n*sizeof(T); i++) {

            // form bitwise xor
            unsigned char bwo = px[i]^py[i];

            // get number of bits from LUT
            result += m_lut[bwo];

        }

        return (double)(8*m_n*sizeof(T)-result);           // FIXME: check for overflow

    }

private:

    using CMercerKernel<T>::m_n;

};


}

#endif
