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

#ifndef NABLA_H
#define NABLA_H

#include "types.h"
#include "lm.h"


class CImageDenoising:public R4R::CLMTVMatrices<R4R::smatf,float> {

public:

    //! Constructor.
    CImageDenoising(size_t width, size_t height):m_width(width),m_height(height){}

    //! \copydoc CLMTVMatrices::ComputeJacobian(R4R::smatf&)
    void ComputeJacobian(R4R::smatf& J);

    //! \copydoc CLMTVMatrices::ComputeGradientOperator(R4R::smatf&)
    void ComputeGradientOperator(R4R::smatf& nabla);

private:

    size_t m_width;
    size_t m_height;

};

#endif // NABLA_H
