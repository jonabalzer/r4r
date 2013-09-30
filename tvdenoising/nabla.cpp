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

#include "nabla.h"

using namespace R4R;

void CImageDenoising::ComputeJacobian(smatf& J) {

    J = smatf(m_height*m_width,m_height*m_width);

    for(size_t i=0; i<m_height*m_width; i++)
        J.Set(i,i,1);

}

void CImageDenoising::ComputeGradientOperator(smatf& nabla) {

    nabla = smatf(2*m_height*m_width,m_height*m_width);





}
