#include "splinecurve.h"

#include <limits>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

CSplineCurve::CSplineCurve():
    m_d(0),
    m_p(0),
    m_n(0),
    m_k(0),
    m_knot(),
    m_cv() {}

CSplineCurve::CSplineCurve(size_t d, size_t p, size_t n):
    m_d(d),
    m_p(p),
    m_n(n),
    m_k(n+2*(p-1)),
    m_knot((float*)_mm_malloc((n+2*(p-1))*sizeof(float),16),CSplineDeallocator()),
    m_cv((float*)_mm_malloc((n*d)*sizeof(float),16),CSplineDeallocator()) {

}

CSplineCurve::~CSplineCurve() {

    // nothing to do, smart pointer will take of themselves

}


void CSplineCurve::MakeClampedUniformKnotVector(float a, float b) {

    assert(m_k);

    float dt = (b - a)/(m_n - 1);

    for(size_t i=0; i<(m_p - 1); i++) {

        m_knot.get()[i] = a;
        m_knot.get()[m_k-i-1] = b;

    }

    for(size_t i=0; i<m_n; i++)
        m_knot.get()[m_p-1+i] = a + i*dt;

}

void CSplineCurve::Print() {

    cout << "Degree: " << m_p << endl;
    cout << "Dimension: " << m_d << endl;
    cout << "Number of CP: " << m_n << endl;
    cout << "Knots: " << endl;

    for(size_t i=0; i<m_k; i++)
        cout << m_knot.get()[i] << endl;

}


int CSplineCurve::GetSpan(float t) {

    // hint
    float dt = (m_knot.get()[m_k-1] - m_knot.get()[0])/(m_n-1);

    int span = (int)(t/dt);  // replace this by the index of knot array?

    if(t>=m_knot.get()[m_p+span-1] && t<m_knot.get()[m_p+span])
        return span;

    return -1;

}


bool CSplineCurve::EvaluateNurbsBasis(int order, const float* knot, float t, float* N) {

    float a0, a1, x, y;
    const float* k0;
    float *t_k, *k_t, *N0;
    const int d = order-1;
    int j, r;

    t_k = (float*)alloca(d<<4);
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
    x = 1.0-sqrt(std::numeric_limits<float>::epsilon());
    if (N[0]>x) {

        if ( N[0] != 1.0 && N[0]<1.0+sqrt(std::numeric_limits<float>::epsilon())) {

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

        if ( N[d] != 1.0 && N[d] < 1.0 + sqrt(std::numeric_limits<float>::epsilon())) {

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

bool CSplineCurve::EvaluateNurbsBasisDerivatives(int order, const float *knot, int der_count, float *N) {

    float dN, c;
    const float *k0, *k1;
    float *a0, *a1, *ptr, **dk;
    int i, j, k, jmax;

    const int d = order - 1;
    const int Nstride = -der_count*order;

    dk = (float**)alloca( (der_count+1) << 3 ); /* << 3 in case pointers are 8 bytes long */
    a0 = (float*)alloca( (order*(2 + ((d+1)>>1))) << 3 ); /* d for a0, d for a1, d*order/2 for dk[]'s and slop to avoid /2 */
    a1 = a0 + order;

    /* initialize reciprocal of knot differences */
    dk[0] = a1 + order;

    for (k = 0; k < der_count; k++) {

        j = d-k;
        k0 = knot++;
        k1 = k0 + j;

        for (i = 0; i < j; i++)
            dk[k][i] = 1.0/(*k1++ - *k0++);

        dk[k+1] = dk[k] + j;

    }

    dk--;

    N += order;

    for ( i=0; i<order; i++) {

        a0[0] = 1.0;

        for (k = 1; k <= der_count; k++) {

            dN = 0.0;
            j = k-i;

            if (j <= 0) {

                dN = (a1[0] = a0[0]*dk[k][i-k])*N[i];
                j = 1;

            }

            jmax = d-i;

            if (jmax < k) {

                while (j <= jmax) {

                    dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
                    j++;

                }

            }
            else {

                /* sum j all the way to j = k */
                while (j < k) {

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
    dN = c = (float)d;
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
