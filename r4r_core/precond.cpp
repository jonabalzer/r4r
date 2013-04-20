/*
 * precond.cpp
 *
 *  Created on: Apr 12, 2012
 *      Author: jbalzer
 */



#include "precond.h"
#include <stdio.h>
#include <iostream>

using namespace std;

namespace R4R {

template<class Matrix,class Vector,class Scalar>
CPreconditioner<Matrix,Vector,Scalar>::CPreconditioner(Matrix& A):
	m_A(A){

}

template class CPreconditioner<CSparseBandedArray<double>,CDenseArray<double>,double>;
template class CPreconditioner<CSparseArray<double>,CDenseArray<double>,double>;
template class CPreconditioner<CDenseArray<double>,CDenseArray<double>,double>;
template class CPreconditioner<CDenseArray<double>,CDenseVector<double>,double>;
template class CPreconditioner<CDenseArray<float>,CDenseVector<float>,float>;
template class CPreconditioner<CSparseArray<double>,CDenseVector<double>,double>;

template<class Matrix,class Vector,class Scalar>
CSSORPreconditioner<Matrix,Vector,Scalar>::CSSORPreconditioner(Matrix& A, Scalar omega, bool lower):
	CPreconditioner<Matrix,Vector,Scalar>(A),
	m_omega(omega),
	m_D(A) {

	if(lower) {

		m_L = CSparseLowerTriangularArray<Scalar>(A);

		A.Transpose();
		m_U = CSparseUpperTriangularArray<Scalar>(A);
		A.Transpose();

	}
	else {

		m_U = CSparseUpperTriangularArray<Scalar>(A);

		A.Transpose();
		m_L = CSparseLowerTriangularArray<Scalar>(A);
		A.Transpose();

	}

	// compute M
	m_L.ScaleDiagonal(1/m_omega);
	m_D.Invert();
	m_D.Scale(m_omega/(2.0-m_omega));
	m_U.ScaleDiagonal(1/m_omega);

}

template<class Matrix,class Vector,class Scalar>
void CSSORPreconditioner<Matrix,Vector,Scalar>::Solve(Vector& x, Vector& y) {

	m_L.Solve(x,y);
	Vector temp(x);
	m_D.Solve(temp,x);
	m_U.Solve(x,temp);

}


template class CSSORPreconditioner<CSparseBandedArray<double>,CDenseArray<double>,double>;
template class CSSORPreconditioner<CSparseArray<double>,CDenseArray<double>,double>;
template class CSSORPreconditioner<CSparseArray<double>,CDenseVector<double>,double>;

template<class Matrix,class Vector,class Scalar>
CJacobiPreconditioner<Matrix,Vector,Scalar>::CJacobiPreconditioner(Matrix& A):
	CPreconditioner<Matrix,Vector,Scalar>(A),
	m_D(A) {

}

template<class Matrix,class Vector,class Scalar>
void CJacobiPreconditioner<Matrix,Vector,Scalar>::Solve(Vector& x, Vector& y) {

	m_D.Solve(x,y);


}


template class CJacobiPreconditioner<CSparseBandedArray<double>,CDenseArray<double>,double>;
template class CJacobiPreconditioner<CSparseArray<double>,CDenseArray<double>,double>;
template class CJacobiPreconditioner<CSparseArray<double>,CDenseVector<double>,double>;



}





