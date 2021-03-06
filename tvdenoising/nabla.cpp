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

#include <vector>

using namespace R4R;
using namespace std;

void CImageDenoising::ComputeJacobian(CCSRMatrix<float,size_t>& J) {

    J = CCSRMatrix<float>(m_height*m_width,m_height*m_width);

    J.Eye();

}

void CImageDenoising::ComputeGradientOperator(CCSRMatrix<float,size_t>& nabla) {

    /* this is only done once, so we can assemble it the slow way
     * with triplets */
    vector<CCSRTriple<float,size_t> > entries;

    size_t i,j;

    // matrices are stored in col-major order
    for(i=0; i<m_height-1; i++) {

        for(j=0; j<m_width-1; j++) {

            size_t row = j*m_height + i;

            // dudx
            entries.push_back(CCSRTriple<float,size_t>(row,row,-1.0));
            entries.push_back(CCSRTriple<float,size_t>(row,row+m_height,1.0));

            // dudy
            entries.push_back(CCSRTriple<float,size_t>(m_width*m_height+row,row,-1.0));
            entries.push_back(CCSRTriple<float,size_t>(m_width*m_height+row,row+1,1.0));

        }

    }

   // bottom of image, dudx exists
   i = m_height - 1;
   for(j=0; j<m_width-1; j++) {

       size_t row = j*m_height + i;

       entries.push_back(CCSRTriple<float,size_t>(row,row,-1.0));
       entries.push_back(CCSRTriple<float,size_t>(row,row+m_height,1.0));

   }

   // right edge of image, dudy exists
   j = m_width - 1;
   for(i=0; i<m_height-1; i++) {

       size_t row = j*m_height + i;

       entries.push_back(CCSRTriple<float,size_t>(m_width*m_height+row,row,-1.0));
       entries.push_back(CCSRTriple<float,size_t>(m_width*m_height+row,row+1,1.0));


   }

   nabla = CCSRMatrix<float>(2*m_height*m_width,m_height*m_width,entries);

}
