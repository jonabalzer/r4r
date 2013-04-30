/*
 * kfilter.cpp
 *
 *  Created on: May 3, 2012
 *      Author: jbalzer
 */

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
