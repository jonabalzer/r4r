//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
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
//////////////////////////////////////////////////////////////////////////////////

#ifndef R4RPRECOND_H_
#define R4RPRECOND_H_


#include "sarray.h"
#include "darray.h"


namespace R4R {

/*! \brief preconditioning of iterative linear solvers
 *
 *
 */
template<class Matrix,typename T>
class CPreconditioner {

public:

	//! Performs preconditioning.
    virtual void Solve(CDenseArray<T>& x, const CDenseArray<T>& y) const { x = y; }

protected:

};

/*! \brief successive over-relaxation preconditioner
 *
 *
 */
template<class Matrix,typename T>
class CSSORPreconditioner: public CPreconditioner<Matrix,T> {

public:

	//! Constructor.
    CSSORPreconditioner(Matrix& A, T omega, bool lower = true);

	//! \copydoc CPreconditioner::Solve(Vector& x, Vector& y)
    void Solve(CDenseArray<T>& x, const CDenseArray<T>& y) const;

protected:

    T m_omega;                                              //!< relaxation parameter
    CSparseLowerTriangularArray<T> m_L;						//!< lower-triangular part of #m_A (or transpose of #m_U)
    CSparseDiagonalArray<T> m_D;							//!< diagonal of #m_A
    CSparseUpperTriangularArray<T> m_U;						//!< upper-triangular part of #m_A (or transpose of #m_L)


};



/*! \brief Jacobi preconditioner
 *
 *
 */
template<class Matrix,typename T>
class CJacobiPreconditioner:public CPreconditioner<Matrix,T> {

public:

	//! Constructor.
	CJacobiPreconditioner(Matrix& A);

	//! \copydoc CPreconditioner::Solve(Vector& x, Vector& y)
    void Solve(CDenseArray<T>& x, const CDenseArray<T>& y) const;

protected:

    CSparseDiagonalArray<T> m_D;						//!< diagonal of #m_A

};



}





#endif /* PRECOND_H_ */
