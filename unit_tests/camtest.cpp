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

#include "camtest.h"

using namespace R4R;

CCamTest::CCamTest(QObject* parent):
  QObject(parent){

}

void CCamTest::init() {

    m_cam = new CPinholeCam<float>(640,480,500,500,319.5,239.5);
    m_view = new CView<float>(*m_cam);

    CRigidMotion<float,3> F(1.2,0,0,0,0,0);
    m_view->SetTransformation(F);

}

void CCamTest::testModelViewProjectionMatrix() {

  matf PF = m_view->ModelViewProjectionMatrix(0.1,100.0);

  float tolerance = 1e-6;

  QVERIFY(tolerance>fabs(1.5625-PF.Get(0,0)));
  QVERIFY(tolerance>fabs(PF.Get(0,1)));
  QVERIFY(tolerance>fabs(-0.0015625-PF.Get(0,2)));
  QVERIFY(tolerance>fabs(1.875-PF.Get(0,3)));

}

void CCamTest::cleanup(){

  delete m_view;
  delete m_cam;

}
