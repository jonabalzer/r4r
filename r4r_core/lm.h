/*
 * lm.h
 *
 *  Created on: May 30, 2012
 *      Author: jbalzer
 */

#ifndef LM_H_
#define LM_H_

#include "iter.h"
#include "types.h"
#include <vector>

namespace R4R {

/*! \brief interface for least-squares problems
 *
 *
 *
 */
template<class Matrix>
class CLeastSquaresProblem {

public:

	//! \brief Standard constructor.
	CLeastSquaresProblem();

	/*! \brief Constructor.
	 *
	 * \param[in] nopts number of data points
	 * \param[in] noparams number of model parameters
	 *
	 */
	CLeastSquaresProblem(size_t nopts, size_t noparams);

	//! \brief Jointly computes the residual vector and Jacobian of the least-squares objective function.
	virtual void ComputeResidualAndJacobian(vec& r, Matrix& J) {};

	//! Computes the residual vector of a weighted least-squares objective function.
	virtual void ComputeResidual(vec& r) {};

	//! Access to the model parameters.
	vec& Get() { return m_model; };

	//! Access to the weights.
	vec& GetWeights() { return m_weights; };

	//! Access to #m_nopts.
	size_t GetNumberOfDataPoints() { return m_nopts; };

	//! Access to #m_noparams.
	size_t GetNumberOfModelParameters() { return m_noparams; };

	/*! \brief Computes scattering of residuals to normalize them for re-weighting and/or outlier detection.
	 *
	 * \details Implementation of the base class assumes univariate data and returns a vector filled with the median
	 * absolute deviation (MAD). Multivariate data require more sophisticated scale-location estimation. In that case,
	 * overload this with implementation of e.g. techniques such as the Donoho-Stahel or Minimum Covariance Determinant
	 * (MCD) estimator.
	 */
	virtual vec ComputeDispersion(const vec& r);

protected:

	size_t m_nopts;								//!< number of data points
	size_t m_noparams;							//!< number of unknown model parameters
	vec m_model;								//!< unknown model parameters
	vec m_weights;								//!< weights

};


/*! \brief Levenberg-Marquardt algorithm to solve nonlinear least-squares problems
 *
 *
 * \details The part for robust regression implements the following weight functions:
 * - bi-square: \f$w_i(r_i)=\chi_{[-1,1]}(1-r_i^2)^2\f$
 * - Huber: \f$w_i(r_i)=\frac{1}{\max(1,|r_i|)}\f$
 *
 */
template<class Matrix>
class CLevenbergMarquardt {

public:

	//! Constructor.
	CLevenbergMarquardt(CLeastSquaresProblem<Matrix>& problem, CIterativeSolver<Matrix,vec,double>& solver, double tau = 1.0);

	//! Triggers execution of Levenberg-Marquardt steps.
	vec Iterate(size_t n, double epsilon1, double epsilon2, bool silent = true);

	//! Starts robust re-weighted Levenberg-Marquardt algorithm.
	vec Iterate(size_t nouter, size_t ninner, double epsilon, bool silentinner, bool silentouter);

protected:

	CLeastSquaresProblem<Matrix>& m_problem;						//!< least-squares problem
	CIterativeSolver<Matrix,vec,double>& m_solver;					//!< linear solver
	double m_tau;													//!< initial damping parameter weight
	double m_lambda;												//!< damping parameter
	std::vector<double> m_residuals;								//!< residuals
	static const double m_params[5];								//!< parameters \f$\rho_1,\rho_2,\beta,\frac{1}{\gamma},\tau,p\f$, cf. [Nielsen1999]

	//! Computes weights based on bi-square function.
	vec BiSquareWeightFunction(const vec& r, vec& w);

	//! Computes weights based on Huber function.
	vec HuberWeightFunction(const vec& r, vec& w);

};

}

#endif /* LM_H_ */
