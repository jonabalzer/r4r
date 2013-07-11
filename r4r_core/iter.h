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

#ifndef R4RITER_H_
#define R4RITER_H_

#include "sarray.h"
#include  "darray.h"
#include "precond.h"

namespace R4R {

/*! \brief iterative solvers for linear equation systems
 *
 * \todo Pre-conditioning in CGLS.
 *
 */
template<class Matrix,class Vector,class Scalar>
class CIterativeSolver {

public:


	/*! \brief Constructor.
	 *
	 * \param[in] M preconditioner
	 * \param[in] n (maximum) number of iterations
	 * \param[in] eps absolute residual at which to terminate
	 * \param[in] silent display on/off
	 *
	 */
	CIterativeSolver(CPreconditioner<Matrix,Vector,Scalar>& M, size_t n, Scalar eps, bool silent = true);

	//! Access to accuracy.
    void SetEps(Scalar eps) { m_eps = eps; }

	//! Access to maximum number of iterations.
    void SetN(size_t n) { m_n = n; }

	/*! \brief Conjugate gradient method.
	 *
	 * \details Implements the conjugate gradient method, cf. [Hestenes1952]. It is crucial that the system matrix \f$A\f$
	 * be symmetric and positive-definite.
	 *
	 * \param[in] A matrix (preserved during operation)
	 * \param[in] b right-hand side
	 * \param[in] x solution vector allocated from the outside
	 * \returns error flag
	 *
	 */
	Scalar CG(Matrix& A, Vector& b, Vector& x);

	//! Conjugate gradient for least-squares method.
	Scalar CGLS(Matrix& A, Vector& b, Vector& x);

	//! Conjugate gradient squared method, cf. [Sonneveld1989].
	Scalar CGS(Matrix& A, Vector& b, Vector& x);

	//! LSQR method, cf. [Paige1982].
    Scalar LSQR(Matrix& A, Vector& b, Vector& x, Scalar lambda = 0);

protected:

    CPreconditioner<Matrix,Vector,Scalar>& m_M;						//!< preconditioner
    size_t m_n;														//!< number of steps
    Scalar m_eps;													//!< absolute accuracy
    bool m_silent;													//!< flag that determines whether messages are displayed

};

}

#endif /* ITER_H_ */
