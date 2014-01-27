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

#include "precond.h"

#include <stdio.h>
#include <iostream>
#include <assert.h>

using namespace std;

namespace R4R {

template class CPreconditioner<CSparseArray<double>,double>;
template class CPreconditioner<CDenseArray<double>,double>;
template class CPreconditioner<CSparseArray<float>,float>;
template class CPreconditioner<CDenseArray<float>,float>;
template class CPreconditioner<CSymmetricCSRMatrix<float,size_t>,float>;
template class CPreconditioner<CSymmetricCSRMatrix<double,size_t>,double>;
template class CPreconditioner<CCSRMatrix<float,size_t>,float>;
template class CPreconditioner<CCSRMatrix<double,size_t>,double>;

template<class Matrix,typename T>
CSSORPreconditioner<Matrix,T>::CSSORPreconditioner(Matrix& A, T omega, bool lower):
    m_D(A),
    m_omega(omega) {

	if(lower) {

        m_L = CSparseLowerTriangularArray<T>(A);

		A.Transpose();
        m_U = CSparseUpperTriangularArray<T>(A);
		A.Transpose();

	}
	else {

        m_U = CSparseUpperTriangularArray<T>(A);

		A.Transpose();
        m_L = CSparseLowerTriangularArray<T>(A);
		A.Transpose();

	}

	// compute M
	m_L.ScaleDiagonal(1/m_omega);
	m_D.Invert();
	m_D.Scale(m_omega/(2.0-m_omega));
	m_U.ScaleDiagonal(1/m_omega);

}

template<class Matrix,typename T>
void CSSORPreconditioner<Matrix,T>::Solve(CDenseArray<T>& x, const CDenseArray<T>& y) const {

	m_L.Solve(x,y);
    CDenseArray<T> temp = x.Clone();
	m_D.Solve(temp,x);
	m_U.Solve(x,temp);

}

template class CSSORPreconditioner<CDenseArray<double>,double>;
template class CSSORPreconditioner<CSparseArray<double>,double>;
template class CSSORPreconditioner<CDenseArray<float>,double>;
template class CSSORPreconditioner<CSparseArray<float>,double>;

template<class Matrix,typename T>
CJacobiPreconditioner<Matrix,T>::CJacobiPreconditioner(Matrix& A):
    m_D(A) {

}

template<class Matrix,typename T>
void CJacobiPreconditioner<Matrix,T>::Solve(CDenseArray<T>& x, const CDenseArray<T>& y) const {

    m_D.Solve(x,y);

}

template class CJacobiPreconditioner<CDenseArray<double>,double>;
template class CJacobiPreconditioner<CSparseArray<double>,double>;
template class CJacobiPreconditioner<CDenseArray<float>,float>;
template class CJacobiPreconditioner<CSparseArray<float>,float>;


} // end of namespace





