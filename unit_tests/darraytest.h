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

#ifndef DARRAYTEST_H
#define DARRAYTEST_H


#include <QtTest/QtTest>

#include "darray.h"

class CDenseArrayTest:public QObject {

  Q_OBJECT

public:

  explicit CDenseArrayTest(QObject* parent = nullptr);

private:

    R4R::CDenseArray<float> m_float_reqc;                //!< array with the same number of rows and columns
    R4R::CDenseArray<float> m_float_rgtc;                //!< more rows than columns
    R4R::CDenseArray<float> m_float_rltc;                //!< more columns than rows
    R4R::CDenseVector<float> m_float_lv;                 //!< larger column vector
    R4R::CDenseVector<float> m_float_sv;                 //!< smaller column vector

private slots:

  void init();

  //! Tests size query.
  void testSize();

  //! Tests matrix-vector multiplication.
  void testMultiplication();

  void cleanup();

};

#endif // DARRAYTEST_H
