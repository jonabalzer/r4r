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

#include "kernelstest.h"

using namespace R4R;
using namespace std;

CKernelsTest::CKernelsTest(QObject* parent):
  QObject(parent) {

}

void CKernelsTest::init() {

    m_n = 1000000;

    m_kernel = new CMercerKernel<float>(m_n);
    m_chi_squared_kernel = new CChiSquaredKernel<float>(m_n);
    m_intersection_kernel = new CIntersectionKernel<float>(m_n);
    m_hellinger_kernel = new CHellingerKernel<float>(m_n);

#ifndef __SSE4_1__
    m_x = new float[m_n];
    m_y = new float[m_n];
#else
    m_x = (float*)_mm_malloc(m_n*sizeof(float),16);
    m_y = (float*)_mm_malloc(m_n*sizeof(float),16);
#endif

    srand(time(NULL));
    for(size_t i=0; i<m_n; i++) {

        m_x[i] = 1.0; // (float)rand()/(float)RAND_MAX+1;
        m_y[i] = 2.0; //(float)rand()/(float)RAND_MAX+1;

    }

}

void CKernelsTest::testIdendityKernel() {

    double parallel = 0;

    QBENCHMARK {

        parallel = m_kernel->Evaluate(m_x,m_y);

    }

    double sequential = 0;

       for(size_t i=0; i<m_n; i++)
            sequential += m_x[i]*m_y[i];

    QCOMPARE(sequential,parallel);

}

void CKernelsTest::testChiSquaredKernel() {

    double parallel = 0;

    QBENCHMARK {

        parallel = m_chi_squared_kernel->Evaluate(m_x,m_y);

    }

    double sequential = 0;

    for(size_t i=0; i<m_n; i++) {

        float num = m_x[i]*m_y[i];

        if(num>0)
            sequential += num/(m_x[i]+m_y[i]);

    }

    QCOMPARE(sequential,parallel);

}


void CKernelsTest::testIntersectionKernel() {

    double parallel = 0;

    QBENCHMARK {

        parallel = m_intersection_kernel->Evaluate(m_x,m_y);

    }

    double sequential = 0;

       for(size_t i=0; i<m_n; i++)
            sequential += min(m_x[i],m_y[i]);

    QCOMPARE(sequential,parallel);

}

void CKernelsTest::testHellingerKernel() {

    double parallel = 0;

    QBENCHMARK {

        parallel = m_hellinger_kernel->Evaluate(m_x,m_y);

    }

    double sequential = 0;

       for(size_t i=0; i<m_n; i++)
            sequential += sqrt(m_x[i]*m_y[i]);

    QCOMPARE(sequential,parallel);

}

void CKernelsTest::cleanup(){

#ifndef __SSE4_1__
    delete [] m_y;
    delete [] m_x;
#else
    _mm_free(m_y);
    _mm_free(m_x);
#endif

    delete m_hellinger_kernel;
    delete m_intersection_kernel;
    delete m_chi_squared_kernel;
    delete m_kernel;

}
