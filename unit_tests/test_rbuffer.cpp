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
#include "rbuffer.h"

using namespace R4R;

class CRingBufferTests:public testing::Test {

protected:

    CRingBufferTests():
        m_float_buffer(5){}

    virtual void SetUp() {

        // [ 3, 5, 1, 0, 0 ]
        m_float_buffer.push_back(3);
        m_float_buffer.push_back(5);
        m_float_buffer.push_back(1);

    }

    CRingBuffer<float> m_float_buffer;              //!< ring buffer of floats

};


TEST_F(CRingBufferTests, BasicFunctionality) {

    EXPECT_EQ(1,m_float_buffer.back());
    EXPECT_EQ(0,m_float_buffer.front());

    // [ 3, 5, 1, 2, 4 ]
    m_float_buffer.push_back(2);
    m_float_buffer.push_back(4);
    EXPECT_EQ(4,m_float_buffer.back());
    EXPECT_EQ(3,m_float_buffer.front());


}
