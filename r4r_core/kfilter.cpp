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

#include "kfilter.h"
#include "factor.h"

namespace R4R {

CKalmanFilter::CKalmanFilter():
	m_x(),
	m_P(),
	m_A(),
	m_At(),
	m_B(),
	m_C(),
	m_Ct(),
	m_Q(),
	m_R() {

}

CKalmanFilter::CKalmanFilter(mat& A, mat& B, mat& C, mat& Q, mat& R):
	m_x(A.NCols()),
	m_P(A.NCols(),A.NCols()),
	m_A(A),
	m_At(A),
	m_B(B),
	m_C(C),
	m_Ct(C),
	m_Q(Q),
	m_R(R) {

	m_At.Transpose();
	m_Ct.Transpose();

}

CKalmanFilter::CKalmanFilter(mat& A, mat& B, mat& C, mat& Q, mat& R, vec& x0):
	m_x(x0),
	m_P(A.NCols(),A.NCols()),
	m_A(A),
	m_At(A),
	m_B(B),
	m_C(C),
	m_Ct(C),
	m_Q(Q),
	m_R(R) {

	m_At.Transpose();
	m_Ct.Transpose();

}


void CKalmanFilter::Predict() {

	m_x = m_A*m_x;

	m_P = m_A*m_P*m_At + m_Q;

}

void CKalmanFilter::Predict(vec& u) {

	m_x = m_A*m_x + m_B*u;

	m_P = m_A*m_P*m_At + m_Q;

}

void CKalmanFilter::Update(vec& y) {

	vec innovation = y - m_C*m_x;

	mat S = m_C*m_P*m_Ct + m_R;

    CMatrixFactorization<double>::InvertSymmetric(S);

	mat K = m_P*m_Ct*S;

	m_x = m_x + K*innovation;

	mat I(m_A.NRows(),m_A.NCols());
	I.Eye();

	m_P = (I - K*m_C)*m_P;

}


}
