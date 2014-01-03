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

#include <Python.h>
#include <numpy/arrayobject.h>
#include <darray.h>
#include <fstream>

#include "io.h"

using namespace std;
using namespace R4R;

static PyMethodDef _iomethods[] = {

    {"read_dense_matrix", read_dense_matrix, METH_VARARGS,"Reads an R4R matrix from file."},
    {NULL, NULL, 0, NULL}

};

PyMODINIT_FUNC init_io(void) {

    (void) Py_InitModule("_io", _iomethods);

    import_array();
}


static PyObject* read_dense_matrix(PyObject* self, PyObject* args) {

    const char* filename;

    // parse tuples
    if (!PyArg_ParseTuple(args,"s",&filename))
        return NULL;

    ifstream in(filename);

    if(!in.is_open()) {

        cerr << "ERROR: File " << filename << " not found..." << endl;
        return NULL;

    }

    // read size
    size_t nrows, ncols;
    in >> nrows;
    in >> ncols;
    int dims[] = { nrows, ncols };

    // read type
    int temp;
    in >> temp;
    ETYPE type = (ETYPE)temp;
    in.get();

    // allocate output
    PyArrayObject* matout;

    switch(type) {

    case ETYPE::B1U:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_BOOL);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(bool));

        break;

    }

    case ETYPE::C1S:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_BYTE);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(char));

        break;

    }

    case ETYPE::C1U:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_UBYTE);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(unsigned char));

        break;

    }

    case ETYPE::S2S:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_SHORT);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(short));

        break;

    }
    case ETYPE::S2U:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_USHORT);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(unsigned short));

        break;

    }

    case ETYPE::I4S:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_INTLTR);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(int));

        break;

    }

    case ETYPE::I4U:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_UINTLTR);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(unsigned int));
        break;

    }

    case ETYPE::F4S:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_FLOAT32);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(float));

        break;

    }

    case ETYPE::L8S:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_LONG);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(long int));

        break;

    }

    case ETYPE::L8U:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_ULONG);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(size_t));

        break;


    }

    case ETYPE::D8S:
    {

        matout = (PyArrayObject*)PyArray_FromDims(2,dims,NPY_DOUBLE);
        in.read((char*)matout->data,(nrows*ncols)*sizeof(double));
        break;

    }
    case ETYPE::D8S3:
    {

        //matout = (PyArrayObject*)PyArray_FromDims(3,dims,NPY_DOUBLE);
        //in.read((char*)matout->data,(nrows*ncols)*sizeof(double));
        break;

    }


    default:
        return NULL;

    }

    return PyArray_Return(matout);

}
