#include "kernels.h"
#include <math.h>
#include <algorithm>

#ifdef __SSE4_1__
#include <xmmintrin.h>
#include <smmintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace R4R {

using namespace std;

template <class T>
double CMercerKernel<T>::Evaluate(T *x, T *y) {

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += (double)x[i]*(double)y[i];

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

template<>
double CMercerKernel<int>::Evaluate(int *x, int *y) {



    return 0;

}

template<class T>
void CMercerKernel<T>::Gradient(T* x, T* y, T* nablax) {

    for(size_t i=0; i<m_n; i++)
        nablax[i] = y[i];

}

template <class T>
CMercerKernel<T>* CMercerKernel<T>::Create(int no, int n) {

    CMercerKernel<T>* kernel;

    switch(no) {

    case IDENTITY:
    {
        kernel = new CMercerKernel<T>(n);
        break;
    }
    case CHISQUARED:
    {
        kernel = new CChiSquaredKernel<T>(n);
        break;
    }
    case INTERSECTION:
    {
        kernel = new CIntersectionKernel<T>(n);
        break;
    }
    case HELLINGER:
    {
        kernel = new CHellingerKernel<T>(n);
        break;
    }
    case SPARSEONESIDED:
    {
        kernel = new COneSidedSparseKernel<T>(n);
        break;
    }
    }

    return kernel;

}

template <class T>
void CMercerKernel<T>::TestKernel(int kn, int n, size_t notests) {

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

    // only for sparse-sided
    int indices[10];
    for(size_t i=0; i<10; i++)
        indices[i] = rand()%n;

    CMercerKernel<float>* kernel;

    if(kn!=SPARSEONESIDED)
        kernel = CMercerKernel<float>::Create(kn,n);
    else
        kernel = CMercerKernel<float>::Create(kn,10);

    double t0, t1;

#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif

    float result;

    if(kn!=SPARSEONESIDED) {

        for(size_t k=0; k<notests; k++)
            result = kernel->Evaluate(x,y);

    }
    else {

        COneSidedSparseKernel<float>* osk = (COneSidedSparseKernel<float>*)(kernel);

        for(size_t k=0; k<notests; k++)
            result = osk->Evaluate(x,indices,y);

    }

#ifdef _OPENMP
    t1 = omp_get_wtime();
    cout << "Time SIMD: " << t1-t0 << endl;
#endif

    cout << result << endl;

    float comparison = 0;

#ifdef _OPENMP
    t0 = omp_get_wtime();
#endif

    switch(kn) {
    {
    case IDENTITY:

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++)
                comparison += x[i]*y[i];

        }

        break;
    }
    case CHISQUARED:
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
    case INTERSECTION:
    {

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++)
                comparison += min(x[i],y[i]);

        }

        break;

    }
    case HELLINGER:
    {

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<n;i++)
                comparison += sqrt(x[i]*y[i]);

        }

        break;

    }
    case SPARSEONESIDED:
    {

        for(size_t k=0; k<notests; k++) {

            comparison = 0;

            for(size_t i=0;i<10;i++)
                comparison += x[i]*y[indices[i]];

        }

    }
    }


#ifdef _OPENMP
    t1 = omp_get_wtime();
    cout << "Time: " << t1-t0 << endl;
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

template <class T>
double COneSidedSparseKernel<T>::Evaluate(T *x, int* indices, T* y) {

    double result = 0;

    for(size_t i=0; i<m_n; i++)
        result += (double)x[indices[i]]*(double)y[i];

    return result;

}

}
