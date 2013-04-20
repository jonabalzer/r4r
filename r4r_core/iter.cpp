/*
 * iter.cpp
 *
 *  Created on: Apr 5, 2012
 *      Author: jbalzer
 */

#include "iter.h"
#include "darray.h"
#include "rutils.h"

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>

namespace R4R {

using namespace std;

template<class Matrix,class Vector,class Scalar>
CIterativeSolver<Matrix,Vector,Scalar>::CIterativeSolver(CPreconditioner<Matrix,Vector,Scalar>& M, size_t n, Scalar eps, bool silent):
	m_M(M),
	m_n(n),
	m_eps(eps),
	m_silent(silent) {

}

template<class Matrix,class Vector,class Scalar>
Scalar CIterativeSolver<Matrix,Vector,Scalar>::CG(Matrix& A, Vector& b, Vector& x) {

	// check dimensions
	if(!(A.NCols()==x.NRows() && x.NRows()==b.NRows() && x.NCols()==b.NCols() && x.NCols()==1)) {

		cout << "ERROR: Check matrix dimensions!" << endl;
		return -1;

	}

	// init
	size_t k = 0;

    Vector r = b - A*x;

	Scalar normr = r.Norm2();

	if(!m_silent)
		cout << "k=" << k << ": " << normr << endl;

	Vector z(r);

	m_M.Solve(z,r);

	Scalar deltan = CDenseArray<Scalar>::InnerProduct(z,r);

	Vector p = z;

	while(k<x.NRows() && k<m_n) {

		Vector q = A*p;

		Scalar alpha = deltan/CDenseArray<Scalar>::InnerProduct(p,q);

		x = x + p*alpha;

		q.Scale(alpha);
		r = r - q;

		normr = r.Norm2();

		k++;

		if(!m_silent)
			cout << "k=" << k << ": " << normr << endl;

		if(normr<m_eps)
			break;

		m_M.Solve(z,r);

		Scalar deltao = deltan;

		deltan = CDenseArray<Scalar>::InnerProduct(z,r);

		Scalar beta = deltan/deltao;

		p.Scale(beta);
		p = z + p;

	}

	return deltan;

}

template<class Matrix,class Vector,class Scalar>
Scalar CIterativeSolver<Matrix,Vector,Scalar>::CGLS(Matrix& A, Vector& b, Vector& x) {

	if(!(A.NCols()==x.NRows() && A.NRows()==b.NRows() && x.NCols()==1)) {

		cout << "ERROR: Check matrix dimensions!" << endl;
		return -1;

	}

	// init
	size_t k = 0;

	// residual of the non-square system
	Vector r = b - A*x;

	Scalar normrt, normr;
	normr = r.Norm2();

	// transpose in place to save computation time
	A.Transpose();
	Vector rnormal = A*r;   // residual of the normal equation
	A.Transpose();


	// preconditioning
	Vector z(rnormal);
	m_M.Solve(z,rnormal);

	// descent direction
	Vector p = z;

	Scalar deltao = CDenseArray<Scalar>::InnerProduct(z,rnormal);  // numerator: z\dot rnormal

	if(!m_silent)
		cout << "k=" << k << ": " << normr << endl;

	//size_t mind = min(A.NRows(),A.NCols());

	//while(k<mind && k<m_n) {
	while(k<m_n) {

		// need that later
		Vector q = A*p;

		// step size (recycle deltao for computation of beta)
		Scalar alpha = deltao/CDenseArray<Scalar>::InnerProduct(q,q);

		// perform descent step
		x = x + p*alpha;

		q.Scale(alpha);
		r = r - q;				// update residual of non-square system

		normrt = r.Norm2();

		k++;

		if(!m_silent)
			cout << "k=" << k << ": " << normr << endl;

		if(fabs(normr-normrt)<m_eps)
			break;

		normr = normrt;

		A.Transpose();
		rnormal = A*r;			// update residual of normal equation
		A.Transpose();

		// apply preconditioner
		m_M.Solve(z,rnormal);
		//z = rnormal;

		// update beta
		Scalar deltan = CDenseArray<Scalar>::InnerProduct(z,rnormal);
		Scalar beta = deltan/deltao;
		deltao = deltan;

		// update direction
		p.Scale(beta);
		p = z + p;

	}

	return deltao;

}

template<class Matrix,class Vector,class Scalar>
Scalar CIterativeSolver<Matrix,Vector,Scalar>::CGS(Matrix& A, Vector& b, Vector& x) {

	// check dimensions
	if(!(A.NCols()==x.NRows() && x.NRows()==b.NRows() && x.NCols()==b.NCols() && x.NCols()==1)) {

		cout << "ERROR: Check matrix dimensions!" << endl;
		return -1;

	}

	// init
	size_t i = 0;

	Vector r = b - A*x;
	Vector rld = r;

	Scalar deltan = r.Norm2();

	if(!m_silent)
		cout << "k=" << i << ": " << deltan << endl;

	Scalar rho1, rho, beta, alpha;
	Vector p, q, u, v, temp;

	while(i<m_n && i<x.NRows() && deltan>m_eps) {

		rho = CDenseArray<Scalar>::InnerProduct(rld,r);

		if(rho==0)
			return -1;

		if (i>0) {

			beta = rho/rho1;

			q.Scale(beta);
			u = r + q;

			temp = u + q;

			Scalar betas = beta*beta;

			p.Scale(betas);

			p = temp + p;

		} else {

			u = r;
			p = u;

		}

		Vector vhat(p);

		m_M.Solve(vhat,p);

		v = A*vhat;

		alpha = rho/CDenseArray<Scalar>::InnerProduct(rld,v);

		v.Scale(alpha);

		q = u - v;

		temp = u + q;

		m_M.Solve(u,temp);

		u.Scale(alpha);

		x = x + u;

		r = r - A*u;

		i++;

		deltan = r.Norm2();

		if(!m_silent)
			cout << "k=" << i << ": " << deltan << endl;

		rho1 = rho;

	}

	return deltan;

}


template<class Matrix,class Vector,class Scalar>
Scalar CIterativeSolver<Matrix,Vector,Scalar>::LSQR(Matrix& A, Vector& b, Vector& x, Scalar lambda) {

	// check dimensions
	if(!(A.NCols()==x.NRows() && x.NRows()==b.NRows() && x.NCols()==b.NCols() && x.NCols()==1)) {

		cout << "ERROR: Check matrix dimensions!" << endl;
		return -1;

	}

	// initialization
	size_t i = 0;

	// temp vector for applying preconditioner
	Vector z = x;

	Vector u = b - A*x;

	Scalar beta = u.Norm2();

	if(beta!=0)
		u.Scale(1/beta);

	m_M.Solve(z,u);

	A.Transpose(); 		// in-place transpose to avoid making expensive copy of A
	Vector v = A*z;
	A.Transpose();

	Scalar alpha = v.Norm2();

	if(alpha!=0)
		v.Scale(1/alpha);

	Vector w = v;
	Scalar phibar = beta;
	Scalar rhobar = alpha;

	Scalar deltan = phibar;

	// check for premature termination
	if(alpha*beta==0)
		return 0;

	if(!m_silent)
		cout << "k=" << i << ": " << deltan << endl;

	while(i<m_n && i<x.NRows()) {

		// 1. Lanczos step
		u.Scale(alpha);

		m_M.Solve(z,v);
		u = A*z - u;

		beta = u.Norm2();

		if(beta!=0)
			u.Scale(1/beta);

		// update v
		v.Scale(beta);

		m_M.Solve(z,u);

		A.Transpose();
		v = A*z - v;
		A.Transpose();

		alpha = v.Norm2();

		if(alpha!=0)
			v.Scale(1/alpha);

		/*Scalar cl, sl, rl;
		R4R::CLinearAlgebra::GivensRotation(rhobar,lambda,cl,sl,rl);
		phibar = cl*phibar;
		Scalar temp = sqrt(rhobar*rhobar + lambda*lambda);*/

		// 2. solve sub least-squares problem by cheap QR decomposition
		Scalar c, s, r;

		CLinearAlgebra::GivensRotation(rhobar,beta,c,s,r);

		Scalar theta = s*alpha;

		rhobar = -c*alpha;

		Scalar phi = c*phibar;

		phibar = s*phibar;

		// 3. update x
		x = x + w*(phi/r);
		w = v - w*(theta/r);

		i++;

		// 4. check termination criteria
		deltan = phibar;

		if(!m_silent)
			cout << "k=" << i << ": " << deltan << endl;

		if(deltan<m_eps)
			break;

	}

	return deltan;

}

template class  CIterativeSolver<CDenseArray<double>,CDenseArray<double>,double>;
template class  CIterativeSolver<CDenseArray<double>,CDenseVector<double>,double>;
template class  CIterativeSolver<CDenseArray<float>,CDenseVector<float>,float>;
template class  CIterativeSolver<CSparseBandedArray<double>,CDenseArray<double>,double>;
template class  CIterativeSolver<CSparseArray<double>,CDenseArray<double>,double>;
template class  CIterativeSolver<CSparseArray<double>,CDenseVector<double>,double>;


}
