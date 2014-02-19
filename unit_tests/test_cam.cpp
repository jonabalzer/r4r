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
#include "cam.h"

using namespace R4R;
using namespace std;

class CViewTests:public testing::Test {

protected:

    CViewTests():
        m_cam(640,480,500,500,319.5,239.5),
        m_view(m_cam) {}



    virtual void SetUp() {

        CRigidMotion<float,3> F(1.2,0,0,0,0,0);
        m_view.SetTransformation(F);


    }

    CPinholeCam<float> m_cam;               //! intrinsic parameters
    CView<float> m_view;                    //! extrinsic parameters

};

TEST_F(CViewTests, OpenGLMatrix) {

    matf PF = m_view.ModelViewProjectionMatrix(0.1,100.0);

    float tolerance = 1e-3;

    EXPECT_GE(tolerance,fabs(1.5625-PF.Get(0,0)));
    EXPECT_GE(tolerance,fabs(PF.Get(0,1)));
    EXPECT_GE(tolerance,fabs(-0.0015625-PF.Get(0,2)));
    EXPECT_GE(tolerance,fabs(1.875-PF.Get(0,3)));

}
