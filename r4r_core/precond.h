/*
 * precond.h
 *
 *  Created on: Apr 12, 2012
 *      Author: jbalzer
 */

#ifndef PRECOND_H_
#define PRECOND_H_


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
