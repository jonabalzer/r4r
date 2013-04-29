#include "splinecurve.h"

#include <limits>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

namespace R4R {

template <class T>
CSplineCurve<T>::CSplineCurve():
    m_d(0),
    m_p(0),
    m_n(0),
    m_k(0),
    m_knot(),
    m_cv() {}

template <class T>
CSplineCurve<T>::CSplineCurve(u_int d, u_int p, u_int n):
    m_d(d),
    m_p(p),
    m_n(n),
    m_k(n+p-1),
    m_knot(n+p-1),
    m_cv(d,n) {         // make sure that the cps are in column-major ordering for fast access!

}

template <class T>
CSplineCurve<T>::CSplineCurve(const CDenseVector<T>& knot, const CDenseArray<T>& cv):
    m_d(cv.NCols()),
    m_p((knot.NElems()-cv.NRows()+1)/2),
    m_n(cv.NRows()),
    m_k(knot.NElems()),
    m_knot(knot),
    m_cv(cv)
{


}

template <class T>
void CSplineCurve<T>::MakeClampedUniformKnotVector(T a, T b) {

    assert(m_k);

    T dt = (b - a)/(m_k - 2*(m_p-1) - 1);
    cout <<"L" << m_knot.NElems() << endl;
    for(size_t i=0; i<(m_p - 1); i++) {

        m_knot(i) = a;
        m_knot(m_k-i-1) = b;

    }

    for(size_t i=0; i<m_k - 2*(m_p-1); i++)
        m_knot(m_p-1+i) = a + i*dt;

}

template <class T>
void CSplineCurve<T>::Print() {

    cout << "Degree: " << m_p << endl;
    cout << "Dimension: " << m_d << endl;
    cout << "Number of CP: " << m_n << endl;
    cout << "Knots: " << endl;

    for(size_t i=0; i<m_k; i++)
        cout << m_knot.Get(i) << endl;

}


template <class T>
int CSplineCurve<T>::GetSpan(T t) {


    if(t<m_knot.Get(0) || t>m_knot.Get(m_knot.NElems()-1))
        return -1;

    // hint
    T dt = (m_knot.Get(m_k-1) - m_knot.Get(0))/(m_n-1);

    int span = (u_int)(t/dt);  // replace this by the index of knot array, binary search with hint!

    return span;

}


template <class T>
CDenseVector<T> CSplineCurve<T>::Evaluate(T t) {

    u_int order = m_p + 1;

    CDenseVector<T> result(m_d);

    int span = GetSpan(t);
    if(span<0)
        return result;

    T* N = new T[order*order];

    EvaluateNurbsBasis(order,m_knot.Data().get()+span,t,N);

    for(size_t i=0; i<order; i++) {

        // get column
        CDenseVector<T> col = m_cv.GetColumn(span+i).Clone();

        // in-place scaling by b(span+i)
        col.Scale(N[i]);

        // in-place add
        result.Add(col);

    }

    return result;

}

template <class T>
CDenseArray<T> CSplineCurve<T>::Tangent(T t) {

    assert(m_p>=1);

    u_int order = m_p + 1;

    // allocate result and get column views
    CDenseArray<T> result(m_d,2);
    int span = GetSpan(t);
    if(span<0)
        return result;

    CDenseVector<T> x = result.GetColumn(0);
    CDenseVector<T> xt = result.GetColumn(1);
    T* N = new T[order*order];

    EvaluateNurbsBasis(order,m_knot.Data().get()+span,t,N);
    EvaluateNurbsBasisDerivatives(order,m_knot.Data().get()+span,1,N);

    for(size_t i=0;i<order; i++) {

        for(size_t j=0; j<order; j++)
            cout << N[i*order+j] << " ";

        cout << endl;

    }

    for(size_t i=0; i<order; i++) {

        // get copy of cv cols
        CDenseVector<T> colx = m_cv.GetColumn(span+i).Clone();
        CDenseVector<T> colxt = m_cv.GetColumn(span+i).Clone();

        // in-place scaling by b(span+i)
        colx.Scale(N[i]);
        colxt.Scale(N[order+i]);

        // in-place add
        x.Add(colx);
        xt.Add(colxt);

    }

    return result;

}


template <class T>
CDenseArray<T> CSplineCurve<T>::Normal(T t) {

    assert(m_p>=2);

    u_int order = m_p + 1;

    // allocate result and get column views
    CDenseArray<T> result(m_d,3);
    int span = GetSpan(t);
    if(span<0)
        return result;

    CDenseVector<T> x = result.GetColumn(0);
    CDenseVector<T> xt = result.GetColumn(1);
    CDenseVector<T> xtt = result.GetColumn(2);
    T* N = new T[order*order];

    EvaluateNurbsBasis(order,m_knot.Data().get()+span,t,N);
    EvaluateNurbsBasisDerivatives(order,m_knot.Data().get()+span,2,N);

    for(size_t i=0; i<order; i++) {

        // get copy of cv cols
        CDenseVector<T> colx = m_cv.GetColumn(span+i).Clone();
        CDenseVector<T> colxt = m_cv.GetColumn(span+i).Clone();
        CDenseVector<T> colxtt = m_cv.GetColumn(span+i).Clone();

        // in-place scaling by b(span+i)
        colx.Scale(N[i]);
        colxt.Scale(N[order+i]);
        colxtt.Scale(N[order*order+i]);

        // in-place add
        x.Add(colx);
        xt.Add(colxt);
        xtt.Add(colxtt);

    }

    return result;

}


template <class T>
bool CSplineCurve<T>::EvaluateNurbsBasis(u_int order, const T* knot, T t, T* N) {

    T a0, a1, x, y;
    const T* k0;
    T *t_k, *k_t, *N0;
    const u_int d = order-1;
    u_int j, r;

    t_k = (T*)alloca(d<<4);
    k_t = t_k + d;

    if (knot[d-1] == knot[d]) {

        memset( N, 0, order*order*sizeof(*N));
        return 1;

    }

    N  += order*order-1;
    N[0] = 1.0;
    knot += d;
    k0 = knot - 1;

    for (j = 0; j < d; j++ ) {

        N0 = N;
        N -= order+1;
        t_k[j] = t - *k0--;
        k_t[j] = *knot++ - t;

        x = 0.0;

        for (r = 0; r <= j; r++) {

            a0 = t_k[j-r];
            a1 = k_t[r];
            y = N0[r]/(a0 + a1);
            N[r] = x + a1*y;
            x = a0*y;

        }

        N[r] = x;

    }

    //   When t is at an end knot, do a check to
    //   get exact values of basis functions.
    //   The problem being that a0*y above can
    //   fail to be one by a bit or two when knot
    //   values are large.
    x = 1.0-sqrt(std::numeric_limits<T>::epsilon());
    if (N[0]>x) {

        if (N[0]!=1.0 && N[0]<1.0+sqrt(std::numeric_limits<T>::epsilon())) {

            r = 1;

            for (j=1; j<=d && r; j++) {

                if (N[j]!=0.0)
                    r = 0;

            }

            if (r)
                N[0] = 1.0;

        }

    }
    else if (N[d]>x) {

        if (N[d]!=1.0 && N[d]<1.0+sqrt(std::numeric_limits<T>::epsilon())) {

            r = 1;

            for (j=0; j<d && r; j++) {

                if(N[j]!=0.0)
                    r = 0;

            }

            if (r)
                N[d] = 1.0;

        }

    }

    return 0;

}

template <class T>
bool CSplineCurve<T>::EvaluateNurbsBasisDerivatives(u_int order, const T* knot, u_int der_count, T* N) {

    T dN, c;
    const T *k0, *k1;
    T *a0, *a1, *ptr, **dk;
    u_int i, j, k, jmax;

    const u_int d = order - 1;
    const int Nstride = -der_count*order;

    dk = (T**)alloca( (der_count+1) << 3 ); /* << 3 in case pointers are 8 bytes long */
    a0 = (T*)alloca( (order*(2 + ((d+1)>>1))) << 3 ); /* d for a0, d for a1, d*order/2 for dk[]'s and slop to avoid /2 */
    a1 = a0 + order;

    /* initialize reciprocal of knot differences */
    dk[0] = a1 + order;

    for (k=0; k<der_count; k++) {

        j = d-k;
        k0 = knot++;
        k1 = k0 + j;

        for (i=0; i<j; i++)
            dk[k][i] = 1.0/(*k1++ - *k0++);

        dk[k+1] = dk[k] + j;

    }

    dk--;

    N += order;

    for (i=0; i<order; i++) {

        a0[0] = 1.0;

        for (k=1; k<=der_count; k++) {

            dN = 0.0;
            j = k-i;

            if (j<=0) {

                dN = (a1[0] = a0[0]*dk[k][i-k])*N[i];
                j = 1;

            }

            jmax = d-i;

            if (jmax<k) {

                while (j<=jmax) {

                    dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
                    j++;

                }

            }
            else {

                /* sum j all the way to j = k */
                while (j<k) {

                    dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
                    j++;

                }

                dN += (a1[k] = -a0[k-1]*dk[k][i])*N[i+k];

            }

            /* d!/(d-k)!*dN = value of k-th derivative */
            N[i] = dN;
            N += order;

            /* a1[] s for next derivative = linear combination
             * of a[]s used to compute this derivative.
             */
            ptr = a0; a0 = a1; a1 = ptr;

        }

        N += Nstride;

    }

    /* apply d!/(d-k)! scaling factor */
    dN = c = (T)d;
    k = der_count;

    while (k--) {

        i = order;

        while (i--)
            *N++ *= c;

        dN -= 1.0;
        c *= dN;

    }

    return 0;

}

template class CSplineCurve<float>;
template class CSplineCurve<double>;

}
