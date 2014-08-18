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

#ifndef RBUFFERTEST_H
#define RBUFFERTEST_H

#include <QtTest/QtTest>

#include "rbuffer.h"

class CRingBufferTest:public QObject {

  Q_OBJECT

public:

  explicit CRingBufferTest(QObject* parent = nullptr);

private:

  R4R::CRingBuffer<float>* m_buffer;

private slots:

  void init();

  void testPushBack();

  void testIterators();

  void cleanup();

};

#endif // RBUFFERTEST_H
