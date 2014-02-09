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


#include <gtest/gtest.h>
#include "darray.h"

using namespace R4R;

class CDenseArrayTests : public ::testing::Test {

protected:

    CDenseArrayTests():
        m_float_reqc(3,3),
        m_float_rgtc(3,2),
        m_float_rltc(2,3),
        m_float_lv(3),
        m_float_sv(2)
        {}

    virtual void SetUp() {

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
         *
         */
        m_float_sv(0) = 9;
        m_float_sv(1) = 3;

    }

    //virtual void TearDown() {}

    CDenseArray<float> m_float_reqc;                //!< array with the same number of rows and columns
    CDenseArray<float> m_float_rgtc;                //!< more rows than columns
    CDenseArray<float> m_float_rltc;                //!< more columns than rows
    CDenseVector<float> m_float_lv;                 //!< larger column vector
    CDenseVector<float> m_float_sv;                 //!< smaller column vector

};

TEST_F(CDenseArrayTests, SizeQuery) {

    EXPECT_EQ(3, m_float_reqc.NRows());
    EXPECT_EQ(3, m_float_reqc.NCols());
    EXPECT_EQ(3, m_float_rgtc.NRows());
    EXPECT_EQ(2, m_float_rgtc.NCols());
    EXPECT_EQ(2, m_float_rltc.NRows());
    EXPECT_EQ(3, m_float_rltc.NCols());

}

TEST_F(CDenseArrayTests, MatrixVectorMultiplication) {

    CDenseVector<float> x = m_float_rgtc*m_float_sv;

    EXPECT_EQ(27, x.Get(0)) << "Some comment" << std::endl;
    EXPECT_EQ(57, x.Get(1));
    EXPECT_EQ(6, x.Get(2));

}
