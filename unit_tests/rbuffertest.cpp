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

#include "rbuffertest.h"


using namespace R4R;

CRingBufferTest::CRingBufferTest(QObject* parent):
  QObject(parent){

}

void CRingBufferTest::init() {

    // [ 0, 0, 0, 0, 0 ]
    m_buffer = new CRingBuffer<float>(5);

    // [ 3, 5, 1, 0, 0 ]
    m_buffer->push_back(3);
    m_buffer->push_back(5);
    m_buffer->push_back(1);

}

void CRingBufferTest::testPushBack() {

    QCOMPARE(m_buffer->back(),1.0);
    QCOMPARE(m_buffer->front(),0.0);

    // does this affect the member variable?
    m_buffer->push_back(2);
    m_buffer->push_back(4);
    QCOMPARE(m_buffer->back(),4.0);
    QCOMPARE(m_buffer->front(),3.0);

}

void CRingBufferTest::testIterators() {

    // the tests seem to be independent from each other, nice QT
    CRingBuffer<float>::iterator it = m_buffer->begin();
    QCOMPARE(*it,0.0);
    ++it;
    QCOMPARE(*it,0.0);
    ++it;
    QCOMPARE(*it,3.0);

}

void CRingBufferTest::cleanup(){

  delete m_buffer;

}
