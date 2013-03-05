/*
 * kfilter.h
 *
 *  Created on: May 3, 2012
 *      Author: jbalzer
 */

#ifndef KFILTER_H_
#define KFILTER_H_


#include "types.h"

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
	vec GetState() const { return m_x; };

	//! Access to covariance.
	mat GetCovariance() const { return m_P; };


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
