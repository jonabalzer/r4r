#include <Python.h>
#include <numpy/arrayobject.h>

#include "bsplines.h"
#include <float.h>

static PyMethodDef _bsplinesmethods[] = {

    { "cox_de_boor", cox_de_boor, METH_VARARGS, "Evaluates basis functions and their derivatives."},
    { "find_knot_span", find_knot_span, METH_VARARGS, "Finds the knot span of a parameter."},
    { "evaluate_curve", evaluate_curve, METH_VARARGS, "Computes point, tangent, and curvature vector of a spline curve."},
    { "unwrap_phase", unwrap_phase, METH_VARARGS, "Phase unwrapping of 1-d signal."},
    { "characteristic_function", characteristic_function, METH_VARARGS, "Tests whether a point is inside or outside of curve."},
    { NULL, NULL, 0, NULL}

};

PyMODINIT_FUNC init_bsplines(void) {

    (void) Py_InitModule("_bsplines", _bsplinesmethods);

    import_array();
}


static PyObject* cox_de_boor(PyObject* self, PyObject* args) {

    PyArrayObject* knots;
    double t;
    int p;

    // parse tuples
    if (!PyArg_ParseTuple(args,"O!id",&PyArray_Type,&knots,&p,&t))
        return NULL;

    if (knots==NULL) return NULL;

    // allocate output
    int dims[] = { p + 1, p + 1 };
    PyArrayObject* matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_DOUBLE);

    // access to data
    double* pknots = (double*)knots->data;

    // find knot span
    int nknots = knots->dimensions[0];
    int span = -1;
    if(t>=pknots[nknots-1-(p-1)]) {
        span = nknots-1-(p-1);
    }
    else if(t<pknots[p-1])
        span = p-1;
    else
        span = _find_knot_span(t,pknots,p-1,nknots-1-(p-1)-1);

    /* Cox-de-Boor recursion
     *
     * The computation is centered around the knot span. We need p-1 intervals to the
     * left and right, so
     * - the knot pointer is offset by -(p-1)
     * - there are (p-1) additional knots at the end point, either open or periodic
     *
     * TODO: Use the basic formula from the main library.
     *
     */
    _cox_de_boor(p+1,pknots+span-(p-1),t,(double*)matout->data);

    // also compute derivatives
    _cox_de_boor_derivatives(p+1,pknots+span-(p-1),p,(double*)matout->data);

    // return array
    PyObject* tupleresult = PyTuple_New(2);
    PyTuple_SetItem(tupleresult, 0, PyArray_Return(matout));
    PyTuple_SetItem(tupleresult, 1, Py_BuildValue("i",span));

    return tupleresult;


}

static PyObject* find_knot_span(PyObject* self, PyObject* args) {

    PyArrayObject* knots;
    double t;
    int p;

    // parse tuples
    if (!PyArg_ParseTuple(args,"O!id",&PyArray_Type,&knots,&p,&t))
        return NULL;

    if (knots==NULL) return NULL;

    // access to data
    double* cknots = (double*)knots->data;

    // find knot span
    int nknots = knots->dimensions[0];
    int span = -1;
    if(t>=cknots[nknots-1-(p-1)]) {
        span = nknots-1-(p-1);
    }
    else if(t<cknots[p-1])
        span = p-1;
    else
        span = _find_knot_span(t,cknots,p-1,nknots-1-(p-1)-1);

    return Py_BuildValue("i",span);

}

static PyObject* evaluate_curve(PyObject* self, PyObject* args) {

    PyArrayObject* knots;
    PyArrayObject* cp;
    double t;
    int p;

    // parse tuples
    if (!PyArg_ParseTuple(args,"O!O!id",&PyArray_Type,&knots,&PyArray_Type,&cp,&p,&t))
        return NULL;

    if (knots==NULL || cp==NULL) return NULL;

    // access to data
    double* pknots = (double*)knots->data;
    double* pcp = (double*)cp->data;

    // find knot span
    int nknots = knots->dimensions[0];
    int span = -1;
    if(t>=pknots[nknots-1-(p-1)]) {
        span = nknots-1-(p-1);
    }
    else if(t<pknots[p-1])
        span = p-1;
    else
        span = _find_knot_span(t,pknots,p-1,nknots-1-(p-1)-1);

    // allocate memory for basis function values and output
    double* N = new double[(p+1)*(p+1)];
    int dims[] = { cp->dimensions[0] };
    PyArrayObject* x = (PyArrayObject*)PyArray_FromDims(1,dims,NPY_DOUBLE);
    PyArrayObject* xt = (PyArrayObject*)PyArray_FromDims(1,dims,NPY_DOUBLE);
    PyArrayObject* xtt = (PyArrayObject*)PyArray_FromDims(1,dims,NPY_DOUBLE);
    double* px = (double*)x->data;
    double* pxt = (double*)xt->data;
    double* pxtt = (double*)xtt->data;

    // evaluate basis function
    _cox_de_boor(p+1,pknots+span-(p-1),t,N);
    _cox_de_boor_derivatives(p+1,pknots+span-(p-1),2,N);

    for(int i=0; i<=p; i++) {

        // get cp index
        int cpi = (span - (p-1) + i)%cp->dimensions[1];

        for(int j=0; j<cp->dimensions[0]; j++) {

            px[j] = px[j] + pcp[j*cp->dimensions[1]+cpi]*N[i];
            pxt[j] = pxt[j] + pcp[j*cp->dimensions[1]+cpi]*N[(p+1)+i];
            pxtt[j] = pxtt[j] + pcp[j*cp->dimensions[1]+cpi]*N[2*(p+1)+i];

        }

    }

    // return array
    PyObject* tupleresult = PyTuple_New(3);
    PyTuple_SetItem(tupleresult, 0, PyArray_Return(x));
    PyTuple_SetItem(tupleresult, 1, PyArray_Return(xt));
    PyTuple_SetItem(tupleresult, 2, PyArray_Return(xtt));

    // clean up
    delete [] N;

    return tupleresult;

}

static PyObject* unwrap_phase(PyObject* self, PyObject* args) {

    PyArrayObject* x;

    if (!PyArg_ParseTuple(args,"O!",&PyArray_Type,&x))
        return NULL;

    if (x==NULL) return NULL;

    // allocate memory for result
    int dims[] = { x->dimensions[0] };
    PyArrayObject* result = (PyArrayObject*)PyArray_FromDims(1,dims,NPY_DOUBLE);

    // first copy into result
    double* px = (double*)x->data;
    double* presult = (double*)result->data;

    for(int i=0; i<dims[0]; i++)
        presult[i] = px[i];

    _unwrap_phase(presult,dims[0]);

    return PyArray_Return(result);

}

static PyObject* characteristic_function(PyObject* self, PyObject* args) {

    PyArrayObject* x;
    PyArrayObject* y;
    PyArrayObject* curve;

    // parse tuples
    if (!PyArg_ParseTuple(args,"O!O!O!",&PyArray_Type,&x,&PyArray_Type,&y,&PyArray_Type,&curve))
        return NULL;

    if (x==NULL || y==NULL || curve==NULL) return NULL;

    // access to data
    double* px = (double*)x->data;
    double* py = (double*)y->data;
    double* pcurve = (double*)curve->data;

    // dims
    int dims[] = { x->dimensions[0], x->dimensions[1] };
    int n = curve->dimensions[0];

    // allocate memory for result
    PyArrayObject* result = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_DOUBLE);
    double* presult = (double*)result->data;
    double* phi = new double[n];

    double eps = 1.0/double(n);

    for(int i=0; i<dims[0]; i++) {

        for(int j=0; j<dims[1]; j++) {

            for(int k=0; k<n; k++) {

                double yk, xk;
                xk = pcurve[2*k] - px[i*dims[1]+j];
                yk = pcurve[2*k+1] - py[i*dims[1]+j];
                phi[k] = atan2(yk,xk);

            }

            _unwrap_phase(phi,n);

            if(fabs(phi[0]-phi[n-1])>eps)
                presult[i*dims[1]+j] = 1.0;
            else
                presult[i*dims[1]+j] = 0;

        }

    }

    delete [] phi;

    return PyArray_Return(result);

}


void _cox_de_boor(int order, double *knot, double t, double* N) {

    double a0, a1, x, y;
    double* k0;
    double *t_k, *k_t, *N0;
    const int d = order - 1;
    int j, r;

    t_k = (double*)alloca(d<<4);
    k_t = t_k + d;

    if (knot[d-1] == knot[d]) {

        memset(N,0,order*order*sizeof(*N));
        return;

    }

    N  += order*order-1;
    N[0] = 1.0;
    knot += d;
    k0 = knot - 1;

    for (j = 0; j<d; j++) {

        N0 = N;
        N -= order+1;
        t_k[j] = t - *k0--;
        k_t[j] = *knot++ - t;

        x = 0.0;

        for (r=0; r<=j; r++) {

            a0 = t_k[j-r];
            a1 = k_t[r];
            y = N0[r]/(a0 + a1);
            N[r] = x + a1*y;
            x = a0*y;

        }

        N[r] = x;

    }

    /* When t is at an end knot, do a check to get exact values of basis functions.
    The problem being that a0*y above can fail to be one by a bit or two when knot
    values are large. */
    x = 1.0 - sqrt(DBL_EPSILON);
    if (N[0]>x) {

        if (N[0]!=1.0 && N[0]<1.0+sqrt(DBL_EPSILON)) {

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

        if (N[d]!=1.0 && N[d]<1.0+sqrt(DBL_EPSILON)) {

            r = 1;

            for (j=0; j<d && r; j++) {

                if(N[j]!=0.0)
                    r = 0;

            }

            if (r)
                N[d] = 1.0;

        }

    }

}

void _cox_de_boor_derivatives(int order, double *knot, int der_count, double* N) {

    double dN, c;
    double *k0, *k1;
    double *a0, *a1, *ptr, **dk;
    int i, j, k, jmax;

    const int d = order - 1;
    const int Nstride = -der_count*order;

    dk = (double**)alloca( (der_count+1) << 3 );                // << 3 in case pointers are 8 bytes long
    a0 = (double*)alloca( (order*(2 + ((d+1)>>1))) << 3 );      // d for a0, d for a1, d*order/2 for dk[]'s and slop to avoid /2
    a1 = a0 + order;

    // initialize reciprocal of knot differences
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

                // sum j all the way to j = k
                while (j<k) {

                    dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
                    j++;

                }

                dN += (a1[k] = -a0[k-1]*dk[k][i])*N[i+k];

            }

            // d!/(d-k)!*dN = value of k-th derivative
            N[i] = dN;
            N += order;

            // a1[] s for next derivative = linear combination of a[]s used to compute this derivative.
            ptr = a0; a0 = a1; a1 = ptr;

        }

        N += Nstride;

    }

    // apply d!/(d-k)! scaling factor
    dN = c = (double)d;
    k = der_count;

    while (k--) {

        i = order;

        while (i--)
            *N++ *= c;

        dN -= 1.0;
        c *= dN;

    }

}

int _find_knot_span(double t, double* knots, int lower, int upper) {

    int middle = lower + (upper - lower)/2;

    if (t<knots[middle])
        return _find_knot_span(t,knots,lower,middle-1);
    else if (t>knots[middle+1])     // if it is in the knot span right to middle, we are done
        return _find_knot_span(t,knots,middle+1,upper);
    else
        return middle;

}

void _unwrap_phase(double* x, size_t n) {

    for(size_t i=1; i<n; i++) {

        double d = x[i] - x[i-1];

        if(d>M_PI) {

            for(size_t j=i; j<n; j++)
                x[j] -= 2*M_PI;

        }
        else if(d<-M_PI) {

            for(size_t j=i; j<n; j++)
                x[j] += 2*M_PI;

        }

    }

}


