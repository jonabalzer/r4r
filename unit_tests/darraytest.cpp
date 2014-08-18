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

#include "darraytest.h"

using namespace R4R;

CDenseArrayTest::CDenseArrayTest(QObject* parent):
  QObject(parent),
  m_float_reqc(3,3),
  m_float_rgtc(3,2),
  m_float_rltc(2,3),
  m_float_lv(3),
  m_float_sv(2){

}

void CDenseArrayTest::init() {

    /* [3,0
     *  6,1
     *  0,2]
     */
    m_float_rgtc(0,0) = 3;
    m_float_rgtc(1,0) = 6;
    m_float_rgtc(1,1) = 1;
    m_float_rgtc(2,1) = 2;

    /* [9,
     *  3]
     */
    m_float_sv(0) = 9;
    m_float_sv(1) = 3;

}

void CDenseArrayTest::testSize() {

    QCOMPARE(m_float_reqc.NRows(),3ul);
    QCOMPARE(m_float_reqc.NCols(),3ul);
    QCOMPARE(m_float_rgtc.NRows(),3ul);
    QCOMPARE(m_float_rgtc.NCols(),2ul);
    QCOMPARE(m_float_rltc.NRows(),2ul);
    QCOMPARE(m_float_rltc.NCols(),3ul);

}

void CDenseArrayTest::testMultiplication() {

    CDenseVector<float> x = m_float_rgtc*m_float_sv;

    QCOMPARE(x.Get(0),27.0);
    QCOMPARE(x.Get(1),57.0);
    QCOMPARE(x.Get(2),6.0);

}

void CDenseArrayTest::cleanup(){


}

