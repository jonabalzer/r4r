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

#include "lm.h"
#include <math.h>
#include "rutils.h"

using namespace std;

namespace R4R {


template <class Matrix,typename T>
CLeastSquaresProblem<Matrix,T>::CLeastSquaresProblem():
	m_nopts(0),
	m_noparams(0),
	m_model(0),
	m_weights(0) {


}

template <class Matrix,typename T>
CLeastSquaresProblem<Matrix,T>::CLeastSquaresProblem(size_t nopts, size_t noparams):
	m_nopts(nopts),
	m_noparams(noparams),
	m_model(noparams),
	m_weights(nopts) {

	m_weights.Ones();

}

template <class Matrix,typename T>
CDenseVector<T> CLeastSquaresProblem<Matrix,T>::ComputeDispersion(CDenseVector<T>& r) {

    T fac = (1.0/1.4826);

    T s = fac/r.MAD(); // FIXME: Is this correct????

    CDenseVector<T> sigma(r.NElems());
	sigma.Ones();
	sigma.Scale(s);

	return sigma;

}

template class CLeastSquaresProblem<mat,double>;
template class CLeastSquaresProblem<smat,double>;
template class CLeastSquaresProblem<smatf,float>;

template <class Matrix,typename T>
const T CLevenbergMarquardt<Matrix,T>::m_params[5] = { 0.25, 0.75 , 2, 1.0/3.0, 3.0 };

template <class Matrix,typename T>
CLevenbergMarquardt<Matrix,T>::CLevenbergMarquardt(CLeastSquaresProblem<Matrix,T>& problem, CIterativeSolver<Matrix,CDenseVector<T>,T>& solver, T tau):
	m_problem(problem),
	m_solver(solver),
	m_tau(tau),
	m_lambda(0),
	m_residuals(0) {

}

template <class Matrix,typename T>
CDenseVector<T> CLevenbergMarquardt<Matrix,T>::Iterate(size_t n, T epsilon1, T epsilon2, bool silent) {

	// access to state
    CDenseVector<T>& x = m_problem.Get();

	// initial residual, Jacobian
    CDenseVector<T> r(m_problem.GetNumberOfDataPoints()+m_problem.GetNumberOfModelParameters());
    Matrix J(m_problem.GetNumberOfDataPoints()+m_problem.GetNumberOfModelParameters(),m_problem.GetNumberOfModelParameters());
	m_problem.ComputeResidualAndJacobian(r,J);

	// initial value for lambda, FIXME: do this depending on trace of J'*J
	m_lambda = m_tau*1;

	for(size_t i=0; i<m_problem.GetNumberOfModelParameters(); i++)
		J(m_problem.GetNumberOfDataPoints()+i,i) = sqrt(m_lambda);

	// residual norm
    T res = r.Norm2();
	m_residuals.push_back(res);

	// gradient norm
	J.Transpose();
    CDenseVector<T> grad = J*r;
	J.Transpose();
    T normgrad = grad.Norm2();

	// init other quantities
    T nu = m_params[2];
	size_t k = 0;

	// in verbose mode, print out initial residual, etc.
	if(!silent) {

		cout.setf(ios::scientific,ios::floatfield);

		cout << "k\t f(x)\t ||grad f||\t ||h||\t lambda" << endl;
		cout << k << "\t" << res << "\t" << normgrad << "\t" << 0.0000 << "\t" << m_lambda << endl;

	}

	while(true) {

		// fill in regularization part
		for(size_t i=0; i<m_problem.GetNumberOfModelParameters(); i++)
			J(m_problem.GetNumberOfDataPoints()+i,i) = sqrt(m_lambda);

		// solve linear system
        CDenseVector<T> step(m_problem.GetNumberOfModelParameters());
        m_solver.CGLS(J,r,step);

        // save old state before advancing
        CDenseVector<T> xold = x.Clone();

        // tentative point, everything is allocated so changing pointer is ok
        x = x - step;

        // compute tentative residual
        CDenseVector<T> rt(m_problem.GetNumberOfDataPoints()+m_problem.GetNumberOfModelParameters());
        Matrix Jt(m_problem.GetNumberOfDataPoints()+m_problem.GetNumberOfModelParameters(),m_problem.GetNumberOfModelParameters());
        m_problem.ComputeResidualAndJacobian(rt,Jt);
		res = rt.Norm2();

		// compute rho
        T normstep, descent, dres, dlres, rho;
		dres = m_residuals.back() - res;
        normstep = sqrt(CDenseVector<T>::InnerProduct(step,step));
        descent = -CDenseVector<T>::InnerProduct(step,grad);
		dlres = 0.5*(normstep*normstep*m_lambda - descent);
		rho = dres/dlres;

		// check step size criterion (lambda going to infinity)
		if(normstep<epsilon2*xold.Norm2())
            break;

		// update state if a descent direction is found
        if(rho>0) {

			// it ok now to store the residual norm
			m_residuals.push_back(res);

            // keep Jacobian and residual, should be ok to move pointers
            r = rt;
            J = Jt;

			// update gradient norm
            J.Transpose();
            grad = J*r;
            J.Transpose();
			normgrad = grad.Norm2();

			// push lambda towards Gauss-Newton step
			nu = m_params[2];
            T factor = max(m_params[3],1-(m_params[2]-1)*pow(2*rho-1,m_params[5]));
            m_lambda *= factor;

			k++;

			// print out current state of optimization
			if(!silent)
                cout << k << "\t" << res << "\t" << normgrad << "\t" << normstep << "\t" << m_lambda << endl;

			// check convergence criteria
			if(normgrad<epsilon1 || k==n )
				break;

        }
		else {

			// push towards gradient descent
			if(m_lambda>0)
				m_lambda *= nu;
			else
				m_lambda = 1e-12;

			nu *= 2;

			// restore state because step was unsuccessful
            x = xold;

        }

	}

	return r;

}

template <class Matrix,typename T>
CDenseVector<T> CLevenbergMarquardt<Matrix,T>::Iterate(size_t nouter, size_t ninner, T epsilon, bool silentouter, bool silentinner) {

    vector<T> residuals;

    CDenseVector<T>& weights = m_problem.GetWeights();
    CDenseVector<T> r(m_problem.GetNumberOfDataPoints());
	m_problem.ComputeResidual(r);

    T res = r.Norm2();
	residuals.push_back(res);

    CDenseVector<T> sigma; // = 1;

	size_t k = 0;

	if(!silentouter) {

		cout.setf(ios::scientific,ios::floatfield);
		cout << "k=" << k << ": " << res << endl;

	}

	while(k<nouter) {

		Iterate(ninner,1e-10,1e-10,silentinner);

		// compute weight function from unweighted residual vector residual
		weights.Ones();
		m_problem.ComputeResidual(r);

        // update weights used in inner iteration
		sigma = BiSquareWeightFunction(r,weights);

		// save norm of residual
		res = r.Norm2();
		residuals.push_back(res);

		k++;

		if(!silentouter)
			cout << "k=" << k << ": " << res << endl;

		double dres = res - residuals[residuals.size()-2];

		if(fabs(dres)<epsilon)
			break;

	}

    // express residual in numbers of variance
	for(size_t i=0;i<r.NElems();i++)
		r(i)=r.Get(i)*sigma.Get(i);


	return r;

}

template <class Matrix,typename T>
CDenseVector<T> CLevenbergMarquardt<Matrix,T>::BiSquareWeightFunction(CDenseVector<T>& r, CDenseVector<T>& w) {

	// estimate standard deviation for normalization of residuals
    CDenseVector<T> sigma = m_problem.ComputeDispersion(r);

	// tuning factor from literature
    T tune = 4.685;

	for(size_t i=0; i<m_problem.GetNumberOfDataPoints(); i++) {

        T x = (r.Get(i)*sigma.Get(i))/tune;

		if(x>=-1 && x<=1)
			w(i) = sqrt((1-x*x)*(1-x*x));
		else
			w(i) = 0;

	}

	return sigma;

}

template <class Matrix,typename T>
CDenseVector<T> CLevenbergMarquardt<Matrix,T>::HuberWeightFunction(CDenseVector<T>& r, CDenseVector<T>& w) {

	// estimate standard deviation for normalization of residuals
    CDenseVector<T> sigma = m_problem.ComputeDispersion(r);

	// tuning factor from literature
    T tune = 1.345;

	for(size_t i=0; i<m_problem.GetNumberOfDataPoints(); i++) {

        T x = (r.Get(i)*sigma.Get(i))/tune;

        w(i) = 1/max(1.0,(double)fabs(x));

	}

	return sigma;

}

template class CLevenbergMarquardt<mat,double>;
template class CLevenbergMarquardt<smat,double>;
template class CLevenbergMarquardt<smatf,float>;

}


