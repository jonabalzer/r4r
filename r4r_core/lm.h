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

#ifndef R4RLM_H_
#define R4RLM_H_

#include <vector>

#include "iter.h"
#include "types.h"


namespace R4R {

/*! \brief interface for least-squares problems
 *
 *
 *
 */
template<class Matrix,typename T>
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
    virtual void ComputeResidualAndJacobian(CDenseVector<T>& r, Matrix& J) = 0;

	//! Computes the residual vector of a weighted least-squares objective function.
    virtual void ComputeResidual(CDenseVector<T>& r) = 0;

	//! Access to the model parameters.
    CDenseVector<T>& Get() { return m_model; }

	//! Access to the weights.
    CDenseVector<T>& GetWeights() { return m_weights; }

	//! Access to #m_nopts.
    size_t GetNumberOfDataPoints() { return m_nopts; }

	//! Access to #m_noparams.
    size_t GetNumberOfModelParameters() { return m_noparams; }

	/*! \brief Computes scattering of residuals to normalize them for re-weighting and/or outlier detection.
	 *
     * \details ComputeDispersion() is actually a misnomer. This function computes the inverse of a robust
     * estimate of the variance of residuals, which are later brought to the bi-unit interval by
     * division with \f$\hat{\sigma}$\f (numerically by multiplication.
     *
	 */
    virtual CDenseVector<T> ComputeDispersion(CDenseVector<T>& r);

protected:

	size_t m_nopts;								//!< number of data points
	size_t m_noparams;							//!< number of unknown model parameters
    CDenseVector<T> m_model;					//!< unknown model parameters
    CDenseVector<T> m_weights;					//!< weights

};


/*! \brief Levenberg-Marquardt algorithm to solve nonlinear least-squares problems
 *
 *
 * \details The part for robust regression implements the following weight functions:
 * - bi-square: \f$w_i(r_i)=\chi_{[-1,1]}(1-r_i^2)^2\f$
 * - Huber: \f$w_i(r_i)=\frac{1}{\max(1,|r_i|)}\f$
 *
 */
template<class Matrix,typename T>
class CLevenbergMarquardt {

public:

	//! Constructor.
    CLevenbergMarquardt(CLeastSquaresProblem<Matrix,T>& problem, CIterativeLinearSolver<Matrix,T>& solver, T tau = 1.0);

	//! Triggers execution of Levenberg-Marquardt steps.
    CDenseVector<T> Iterate(size_t n, T epsilon1, T epsilon2, bool silent = true);

	//! Starts robust re-weighted Levenberg-Marquardt algorithm.
    CDenseVector<T> Iterate(size_t nouter, size_t ninner, T epsilon, bool silentinner, bool silentouter);

protected:

    CLeastSquaresProblem<Matrix,T>& m_problem;						//!< least-squares problem
    CIterativeLinearSolver<Matrix,T>& m_solver;                     //!< linear solver
    T m_tau;        												//!< initial damping parameter weight
    T m_lambda;             										//!< damping parameter
    std::vector<T> m_residuals;     								//!< residuals
    static const T m_params[5];                                     //!< parameters \f$\rho_1,\rho_2,\beta,\frac{1}{\gamma},\tau,p\f$, cf. [Nielsen1999]

	//! Computes weights based on bi-square function.
    CDenseVector<T> BiSquareWeightFunction(CDenseVector<T>& r, CDenseVector<T>& w);

	//! Computes weights based on Huber function.
    CDenseVector<T> HuberWeightFunction(CDenseVector<T>& r, CDenseVector<T>& w);

};

/*! \brief Split Bregman method
 *
 */
template<class Matrix,typename T>
class CSplitBregman {

    //! Constructor.
    CSplitBregman(const Matrix& A, const Matrix& Phi, const CDenseArray<T>& f, CDenseArray<T>& u, const CIterativeLinearSolver<Matrix,T>& solver, T mu, T lambda, double eps);

    //! Deleted standard constructor.
    CSplitBregman() = delete;

    /*! \brief Iterate.
     *
     * \param[in] n maximum number of iterations
     *
     */
    void Iterate(size_t n);

    //! Access to residual.
    std::vector<double>& GetTotalError() { return m_total_error; }

    //! Access to constraint violation.
    std::vector<double>& GetConstraintViolation() { return m_constraint_violation; }

private:

    const Matrix& m_K;                                      //! stack of two linear operators
    const Matrix& m_nabla;                                  //! gradient operator
    const CDenseArray<T>& m_f;                              //! force vector
    const CDenseArray<T>& m_u;                              //! solution vector
    const CIterativeLinearSolver<Matrix,T>& m_solver;       //! linear solver
    T m_lambda;                                             //! \f$\lambda\f$
    T m_mu;                                                 //! \f$\mu\f$
    double m_eps;                                           //! tolerance
    std::vector<double> m_total_error;                      //! total error over time
    std::vector<double> m_constraint_violation;             //! constraint violation over time

};

/*! interface for function objects needed in TV-regularized LM method
 *
 */
template <class Matrix,typename T>
class CLMTVMatrices {

public:

    //! Computes residual and Jacobian.
    virtual void ComputeJacobian(Matrix& J) = 0;

    //! Computes representation of gradient operator.
    virtual void ComputeGradientOperator(Matrix& nabla) = 0;

};

}

#endif /* LM_H_ */
