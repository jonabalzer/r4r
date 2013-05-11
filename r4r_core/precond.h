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

#ifndef R4RPRECOND_H_
#define R4RPRECOND_H_


#include "sarray.h"
#include "darray.h"


namespace R4R {

/*! \brief preconditioning of iterative linear solvers
 *
 *
 */
template<class Matrix,class Vector,class Scalar>
class CPreconditioner {

public:

	//! Constructor.
	CPreconditioner(Matrix& A);

	//! Performs preconditioning.
	virtual void Solve(Vector& x, Vector& y) { x = y; };

protected:

	Matrix& m_A;						//!< input matrix, M is member of inherited classes


};

/*! \brief successive over-relaxation preconditioner
 *
 *
 */
template<class Matrix,class Vector,class Scalar>
class CSSORPreconditioner: public CPreconditioner<Matrix,Vector,Scalar> {

public:

	//! Constructor.
	CSSORPreconditioner(Matrix& A, Scalar omega, bool lower = true);

	//! \copydoc CPreconditioner::Solve(Vector& x, Vector& y)
	void Solve(Vector& x, Vector& y);

protected:

	Scalar m_omega;													//!< relaxation parameter
	CSparseLowerTriangularArray<Scalar> m_L;						//!< lower-triangular part of #m_A (or transpose of #m_U)
	CSparseDiagonalArray<Scalar> m_D;								//!< diagonal of #m_A
	CSparseUpperTriangularArray<Scalar> m_U;						//!< upper-triangular part of #m_A (or transpose of #m_L)

	using CPreconditioner<Matrix,Vector,Scalar>::m_A;

};



/*! \brief Jacobi preconditioner
 *
 *
 */
template<class Matrix,class Vector,class Scalar>
class CJacobiPreconditioner:public CPreconditioner<Matrix,Vector,Scalar> {

public:

	//! Constructor.
	CJacobiPreconditioner(Matrix& A);

	//! \copydoc CPreconditioner::Solve(Vector& x, Vector& y)
	void Solve(Vector& x, Vector& y);

protected:

	CSparseDiagonalArray<Scalar> m_D;								//!< diagonal of #m_A

	using CPreconditioner<Matrix,Vector,Scalar>::m_A;

};



}





#endif /* PRECOND_H_ */
