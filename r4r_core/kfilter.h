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

#ifndef R4RKFILTER_H_
#define R4RKFILTER_H_

#include "darray.h"

namespace R4R {

/*! \brief Kalman filter
 *
 *	\details No feed-forward \f$D\f$ implemented, yet.
 *
 */
class CKalmanFilter {

public:

	//! Standard constructor.
	CKalmanFilter();

	//! Constructor.
	CKalmanFilter(mat& A, mat& B, mat& C, mat& Q, mat& R);

	//! Constructor.
	CKalmanFilter(mat& A, mat& B, mat& C,  mat& Q, mat& R, vec& x0);

	//!  Predicts new state.
	void Predict(vec& u);

	//!  Predicts new state without control input.
	void Predict();

	//! Updates state given a measurement.
	void Update(vec& y);

	//! Access to state.
    vec GetState() const { return m_x; }

	//! Access to covariance.
    mat GetCovariance() const { return m_P; }


protected:

	vec m_x;						//!< state
	mat m_P;						//!< state covariance matrix

	mat m_A;						//!< system matrix
	mat m_At;						//!< transpose of system matrix
	mat m_B;						//!< input matrix
	mat m_C;						//!< measurement matrix
	mat m_Ct;						//!< transpose of measurement matrix

	mat m_Q;						//!< system noise covariance matrix
	mat m_R;						//!< measurement noise covariance matrix, TODO: estimate online by ALS

};

}


#endif /* KFILTER_H_ */
