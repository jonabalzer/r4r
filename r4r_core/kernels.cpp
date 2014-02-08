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

#include <math.h>
#include <algorithm>

#ifdef __SSE4_1__
#include <xmmintrin.h>
#include <smmintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kernels.h"
#include "darray.h"


using namespace std;



namespace R4R {

template <class T>
double CMercerKernel<T>::Evaluate(T *x, T *y) {

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += (double)x[i]*(double)y[i];

    return result;

}

template <>
double CMercerKernel<rgb>::Evaluate(rgb* x, rgb* y) {

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += InnerProduct(x[i],y[i]);

    return result;

}

template <>
double CMercerKernel<vec3>::Evaluate(vec3* x, vec3* y) {

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += InnerProduct(x[i],y[i]);

    return result;

}


template<>
double CMercerKernel<float>::Evaluate(float* x, float* y) {

#ifndef __SSE4_1__

    float result = 0;

    for(size_t i=0; i<m_n; i++)
        result += x[i]*y[i];

    return (double)result;

#else
    __m128* px = (__m128*)x;
    __m128* py = (__m128*)y;

    float zero = 0;
    __m128 sum = _mm_load1_ps(&zero);

    const int mask = 241;       // 4 MSB mask input, 4 LSB mask output

    for(size_t i=0; i<m_offset/4; i++) {

        __m128 temp = _mm_dp_ps(px[i],py[i],mask);
        sum = _mm_add_ss(sum,temp);                         // accumulate result in first register

    }

    float result[4] = {0,0,0,0};
    _mm_storeu_ps(result,sum);

    // add offset
    for(size_t i=m_offset; i<m_n; i++)
        result[0] += x[i]*y[i];

    return (double)result[0];
#endif

}

template<class T>
void CMercerKernel<T>::Gradient(T* x, T* y, T* nablax) {

    for(size_t i=0; i<m_n; i++)
        nablax[i] = y[i];

}

template <class T>
CMercerKernel<T>* CMercerKernel<T>::Create(KERNEL no, int n) {

    CMercerKernel<T>* kernel;

    switch(no) {

    case KERNEL::IDENTITY:
    {
        kernel = new CMercerKernel<T>(n);
        break;
    }
    case KERNEL::CHISQUARED:
    {
        kernel = new CChiSquaredKernel<T>(n);
        break;
    }
    case KERNEL::INTERSECTION:
    {
        kernel = new CIntersectionKernel<T>(n);
        break;
    }
    case KERNEL::HELLINGER:
    {
        kernel = new CHellingerKernel<T>(n);
        break;
    }

    }

    return kernel;

}

template <>
void CMercerKernel<float>::TestKernel(KERNEL no, int n, size_t notests) {

    srand(time(NULL));

    cout << "Size (mxn): " << notests << " " << n << endl;

#ifndef __SSE4_1__
    float* x = new float[n];
    float* y = new float[n];
#else
    float* x = (float*)_mm_malloc(n*sizeof(float),16);
    float* y = (float*)_mm_malloc(n*sizeof(float),16);
#endif

    for(size_t i=0; i<n; i++) {

        x[i]=(float)rand()/(float)RAND_MAX+1;
        y[i]=(float)rand()/(float)RAND_MAX+1;

    }

    CMercerKernel<float>* kernel = CMercerKernel<float>::Create(no,n);

#ifdef _OPENMP
    double t0, t1;
    t0 = omp_get_wtime();
#endif

    float result;

    for(size_t k=0; k<notests; k++)
        result = kernel->Evaluate(x,y);

#ifdef _OPENMP
    t1 = omp_get_wtime();
    cout << "Parallel evaluation time: " << t1-t0 << " s" << endl;
#endif

    cout << result << endl;

    float comparison = 0;

#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif

    switch(no) {
    {
    case KERNEL::IDENTITY:

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++)
                comparison += x[i]*y[i];

        }

        break;
    }
    case KERNEL::CHISQUARED:
    {

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++) {

                float num = x[i]*y[i];

                if(num>0)
                    comparison += num/(x[i]+y[i]);


            }
        }

        break;

    }
    case KERNEL::INTERSECTION:
    {

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++)
                comparison += min(x[i],y[i]);

        }

        break;

    }
    case KERNEL::HELLINGER:
    {

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++)
                comparison += sqrt(x[i]*y[i]);

        }

        break;

    }

    }


#ifdef _OPENMP
    t1 = omp_get_wtime();
    cout << "Sequential evaluation time: " << t1-t0 << " s" << endl;
#endif

    cout << comparison << endl;

#ifndef __SSE4_1__
    delete [] x;
    delete [] y;
#else
    _mm_free(y);
    _mm_free(x);
#endif

    delete kernel;

}

template class CMercerKernel<float>;
template class CMercerKernel<double>;
template class CMercerKernel<bool>;
template class CMercerKernel<int>;
template class CMercerKernel<size_t>;
template class CMercerKernel<unsigned char>;

template<class T>
double CChiSquaredKernel<T>::Evaluate(T* x, T* y) {

    double result, num;
    double xi, yi;

    for(size_t i=0; i<m_n; i++) {

        xi = (double)x[i];
        yi = (double)y[i];

        num = xi*yi;

        if(num>0)               // this implies that x+y!=0 if x,y>0
            result += num/(xi+yi);

    }

    return result;

}

template <>
double CChiSquaredKernel<float>::Evaluate(float* x, float* y) {

#ifndef __SSE4_1__

    double result, num;
    double xi, yi;

    for(size_t i=0; i<m_n; i++) {

        xi = (double)x[i];
        yi = (double)y[i];

        num = xi*yi;

        if(num>0)               // this implies that x+y!=0 if x,y>0
            result += num/(xi+yi);

    }

    return result;

#else

    __m128* px = (__m128*)x;
    __m128* py = (__m128*)y;


    __m128 sum = _mm_set1_ps(0.0f);
    __m128 mzero = _mm_set1_ps(0.0f);

    for(size_t i=0; i<m_offset/4; i++) {

        __m128 num = _mm_mul_ps(px[i],py[i]);
        __m128 denom = _mm_add_ps(px[i],py[i]);
        __m128 invdenom = _mm_rcp_ps(denom);

        // find nonzeros in numerator
        __m128 nans = _mm_cmpeq_ps(num,mzero);
        __m128 factors = _mm_blendv_ps(invdenom,mzero,nans);

        // compute product
        __m128 temp = _mm_mul_ps(num,factors);

        // add
        sum = _mm_add_ps(sum,temp);

    }

    float result[4] = {0,0,0,0};
    _mm_storeu_ps(result,sum);

    float fresult  = result[0] + result[1] + result[2] + result[3];

    // add offset
    float num;

    for(size_t i=m_offset; i<m_n; i++)  {

        num = x[i]*y[i];

        if(num>0)                      // this implies that x+y!=0 if x,y>0
            fresult += num/(x[i]+y[i]);

    }

    return (double)fresult;

#endif

}

template class CChiSquaredKernel<float>;
template class CChiSquaredKernel<double>;
template class CChiSquaredKernel<int>;
template class CChiSquaredKernel<size_t>;
template class CChiSquaredKernel<bool>;
template class CChiSquaredKernel<unsigned char>;

template <class T>
double CIntersectionKernel<T>::Evaluate(T* x, T* y) {

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += (double)min<T>(x[i],y[i]);

    return result;

}

template <>
double CIntersectionKernel<float>::Evaluate(float* x, float* y) {

#ifndef __SSE4_1__

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += (double)min<float>(x[i],y[i]);

    return result;

#else

    __m128* px = (__m128*)x;
    __m128* py = (__m128*)y;

    float zero = 0;
    __m128 sum = _mm_load1_ps(&zero);

    const int mask = 255;

    for(size_t i=0; i<m_offset/4; i++) {

        __m128 temp = _mm_min_ps(px[i],py[i]);
        sum = _mm_add_ps(sum,temp);

    }

    float result[4] = {0,0,0,0};
    _mm_storeu_ps(result,sum);

    float fresult  = result[0] + result[1] + result[2] + result[3];

    // add offset
    for(size_t i=m_offset; i<m_n; i++)
        fresult += min<float>(x[i],y[i]);

    return (double)fresult;

#endif

}

template class CIntersectionKernel<float>;
template class CIntersectionKernel<double>;
template class CIntersectionKernel<bool>;
template class CIntersectionKernel<int>;
template class CIntersectionKernel<size_t>;
template class CIntersectionKernel<unsigned char>;

template <class T>
double CHellingerKernel<T>::Evaluate(T* x, T* y) {

    double result = 0;
    double xi, yi;

    for(size_t i=0; i<m_n; i++) {

        xi = (double)x[i];
        yi = (double)y[i];

        result += sqrt(xi*yi);

    }

    return result;

}

template <>
double CHellingerKernel<rgb>::Evaluate(rgb* x, rgb* y) {

    double result = 0;
    double xi, yi;

    for(size_t i=0; i<m_n; i++) {

        rgb xy = x[i]*y[i];

        for(size_t j=0; j<3; j++)
            result += sqrt(xy.Get(j));

    }

    return result;

}

template <>
double CHellingerKernel<vec3>::Evaluate(vec3* x, vec3* y) {

    double result = 0;
    double xi, yi;

    for(size_t i=0; i<m_n; i++) {

        vec3 xy = x[i]*y[i];

        for(size_t j=0; j<3; j++)
            result += sqrt(xy.Get(j));

    }

    return result;

}
template <>
double CHellingerKernel<float>::Evaluate(float* x, float* y) {


#ifndef __SSE4_1__

    double result = 0;
    double xi, yi;

    for(size_t i=0; i<m_n; i++) {

        xi = (double)x[i];
        yi = (double)y[i];

        result += sqrt(xi*yi);

    }

    return result;

#else

    __m128* px = (__m128*)x;
    __m128* py = (__m128*)y;

    float zero = 0;
    __m128 sum = _mm_load1_ps(&zero);

    for(int i=0; i<m_offset/4; i++) {

        __m128 temp = _mm_mul_ps(px[i],py[i]);
        temp = _mm_sqrt_ps(temp);
        sum = _mm_add_ps(sum,temp);

    }

    float result[4] = {0,0,0,0};
    _mm_storeu_ps(result,sum);

    float fresult  = result[0] + result[1] + result[2] + result[3];

    // add offset
    for(size_t i=m_offset; i<m_n; i++)
        fresult += sqrt(x[i]*y[i]);

    return (double)fresult;

#endif

}

template class CHellingerKernel<float>;
template class CHellingerKernel<double>;
template class CHellingerKernel<int>;
template class CHellingerKernel<bool>;
template class CHellingerKernel<size_t>;
template class CHellingerKernel<unsigned char>;

template <class T>
unsigned int CBitwiseHammingKernel<T>::m_lut[256] = {   0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3,
                                                        3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4,
                                                        3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2,
                                                        2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5,
                                                        3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5,
                                                        5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3,
                                                        2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
                                                        4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                                        3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4,
                                                        4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6,
                                                        5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5,
                                                        5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8, };

template class CBitwiseHammingKernel<unsigned char>;
//template class CBitwiseHammingKernel<float>;
template class CHammingKernel<int>;

template <class T>
double CRBFKernel<T>::Evaluate(T *x, T *y) {

    double xx = CMercerKernel<T>::Evaluate(x,x);
    double yy = CMercerKernel<T>::Evaluate(y,y);
    double xy = CMercerKernel<T>::Evaluate(x,y);

    return exp(-(xx-2*xy+xx)/(2*m_sigma*m_sigma));

}

template class CRBFKernel<float>;
template class CRBFKernel<double>;
template class CRBFKernel<bool>;
template class CRBFKernel<size_t>;
template class CRBFKernel<int>;
template class CRBFKernel<unsigned char>;

template <class T>
double CPolynomialKernel<T>::Evaluate(T*x , T*y) {

    double xy = CMercerKernel<T>::Evaluate(x,y);

    return pow(xy+m_c,m_d);

}

template class CPolynomialKernel<float>;
template class CPolynomialKernel<double>;
template class CPolynomialKernel<int>;
template class CPolynomialKernel<bool>;
template class CPolynomialKernel<size_t>;
template class CPolynomialKernel<unsigned char>;

}
