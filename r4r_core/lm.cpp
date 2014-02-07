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

#include <algorithm>
#include <math.h>
#include <limits>
#include <assert.h>

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

    /* this is ok, actually we are computing the inverse of the covariance
     * for normalization of the residuals later, which then only requires multiplication*/
    T fac = (1.0/1.4826);

    T s = fac/r.MAD();

    CDenseVector<T> sigma(r.NElems());
	sigma.Ones();
	sigma.Scale(s);

	return sigma;

}

template class CLeastSquaresProblem<mat,double>;
template class CLeastSquaresProblem<smat,double>;
template class CLeastSquaresProblem<smatf,float>;
template class CLeastSquaresProblem<CCSRMatrix<float>,float>;

template <class Matrix,typename T>
const T CLevenbergMarquardt<Matrix,T>::m_params[5] = { 0.25, 0.75 , 2, 1.0/3.0, 3.0 };

template <class Matrix,typename T>
CLevenbergMarquardt<Matrix,T>::CLevenbergMarquardt(CLeastSquaresProblem<Matrix,T>& problem, CIterativeLinearSolver<Matrix,T>& solver, T tau):
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
    Matrix J(m_problem.GetNumberOfDataPoints(),m_problem.GetNumberOfModelParameters());
	m_problem.ComputeResidualAndJacobian(r,J);

    // initial value for lambda, TODO: do this depending on trace of J'*J
	m_lambda = m_tau*1;

    // residual norm, TODO: why are we not counting the value of the step size?
    T res = r.Norm2();
    m_residuals.push_back(res);

    // init regularization matrix
    Matrix I(m_problem.GetNumberOfModelParameters(),m_problem.GetNumberOfModelParameters());
    I.Eye();
    I.Scale(sqrt(m_lambda));

    // append it to J
    J.Concatenate(I,0);

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

		// solve linear system
        CDenseVector<T> step(m_problem.GetNumberOfModelParameters());
        m_solver.Iterate(J,r,step);

        // save old state before advancing
        CDenseVector<T> xold = x.Clone();

        // tentative point, everything is allocated so changing pointer is ok
        x = x - step;

        // compute tentative residual and Jacobian
        CDenseVector<T> rt(m_problem.GetNumberOfDataPoints()+m_problem.GetNumberOfModelParameters());
        Matrix Jt(m_problem.GetNumberOfDataPoints(),m_problem.GetNumberOfModelParameters());
        m_problem.ComputeResidualAndJacobian(rt,Jt);

        // residual norm, the regularzing part of J is not needed
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

            // add regularizer back
            I.Eye();
            I.Scale(sqrt(m_lambda));
            J.Concatenate(I,0);

            // update gradient norm, this contains step size parameter (but maybe it should not?)
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
            if(normgrad<epsilon1 || k==n)
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

            // show how lambda develops
            if(!silent)
                cout << "Adjusting lambda to " << m_lambda << "." << endl;

        }

        if(std::isinf(m_lambda))
            break;

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

        // run LM method
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
        r(i) = r.Get(i)*sigma.Get(i);


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
template class CLevenbergMarquardt<CCSRMatrix<float>,float>;

template<class Matrix,typename T>
CSplitBregman<Matrix,T>::CSplitBregman(const Matrix& A, const Matrix& nabla, const CDenseArray<T>& f, CDenseArray<T>& u, const DIM& dim, const CIterativeLinearSolver<Matrix,T>& solver, T mu, T lambda, double eps):
    m_K(A.Clone()),
    m_nabla(nabla),
    m_f(f),
    m_u(u),
    m_nablau(),
    m_b(nabla.NRows(),f.NCols()),
    m_d(nabla.NRows(),f.NCols()),
    m_dim_grad(dim),
    m_solver(solver),
    m_mu(mu),
    m_lambda(lambda),
    m_eps(eps),
    m_k(0),
    m_total_error(),
    m_constraint_violation() {

    assert(lambda>0);

    // scale
    m_K.Scale(-mu/lambda);

    // stack on top of m_Phi
    m_K.Concatenate(nabla,0);

    // scale again
    m_K.Scale(-lambda);

    // nabla u and b
    m_nablau = m_nabla*m_u;

    // store first constraint violation: \|\nabla u-d\|=\|\nabla u\|
    m_constraint_violation.push_back(m_nablau.Norm2());

}

template<class Matrix,typename T>
void CSplitBregman<Matrix,T>::Iterate(size_t n) {

    // main loop
    do {

        // create rhs
        CDenseArray<T> rhs(m_K.NRows(),m_f.NCols());

        // fill in m_f and b
        for(size_t i=0; i<m_K.NRows(); i++) {

            for(size_t j=0; j<m_f.NCols(); j++) {

                if(i<m_f.NRows())
                    rhs.Set(i,j,m_f.Get(i,j)*m_mu);
                else
                    rhs.Set(i,j,(m_b.Get(i-m_f.NRows(),j)-m_d.Get(i-m_f.NRows(),j))*m_lambda);

            }

        }

        // get half residual of data term and constraint violation from linear solver
        vector<double> rt = m_solver.Iterate(m_K,rhs,m_u);

        // if k=0, we need compute the initial residual
        if(m_k==0) {

            double rtotal = rt.front()*rt.front(); //mu|Au-f|^2+\lambda\|nabla u-d\|^2
            double cv = m_constraint_violation.back()*m_constraint_violation.back()*m_lambda*m_lambda;
            double temp = 0.5*(rtotal - cv)/m_mu + m_nablau.Norm1();   // remove one mu
            m_total_error.push_back(temp);

            cout << "k" << "\t\t" << "Total error" << "\t\t" << "Feasability" << endl;
            cout << m_k << "\t\t" << m_total_error.back() << "\t\t" << m_constraint_violation.back() << endl;

        }

        // update gradient
        m_nablau = m_nabla*m_u;

        // shrinkage
        m_d = m_b + m_nablau;
        //m_d.Shrink(1/m_lambda);
        this->Shrink();

        // constraint violation
        CDenseArray<T> rphi = m_nablau - m_d;

        // Bregman update
        m_b = m_b + rphi;

        // compute constraint violation
        m_constraint_violation.push_back(rphi.Norm2());

        // total error
        double rtotal = rt.back()*rt.back(); //mu|Au-f|^2+\lambda\|nabla u-d\|^2
        double cv = m_constraint_violation.back()*m_constraint_violation.back()*m_lambda*m_lambda;
        double temp = 0.5*(rtotal - cv)/m_mu + m_nablau.Norm1();   // remove one mu
        m_total_error.push_back(temp);

        // increment
        m_k++;

        // plot errors
        cout << m_k << "\t\t" << m_total_error.back() << "\t\t" << m_constraint_violation.back() << endl;

    } while(m_k<n && fabs(m_total_error.back()-m_total_error.at(m_total_error.size()-2))>m_eps);

}

template<class Matrix,typename T>
void CSplitBregman<Matrix,T>::Shrink() {

    u_int dim = u_int(m_dim_grad);

    assert(m_d.NRows()%dim==0);

    size_t npts = m_d.NRows()/dim;

    T li = 1/m_lambda;

    for(size_t j=0; j<m_d.NCols(); j++) {

        for(size_t i=0; i<npts; i++) {

            T norm = 0;

            for(u_int k=0; k<dim; k++)
                norm += m_d.Get(i+k*npts,j)*m_d.Get(i+k*npts,j);

            norm = sqrt(norm);

            T factor;
            if(norm<=li)
                factor = 0;
            else
                factor = (norm - li)/norm;

            for(u_int k=0; k<dim; k++)
                m_d(i+k*npts,j) *= factor;

        }

    }

}



template class CSplitBregman<smatf,float>;
template class CSplitBregman<CCSRMatrix<float,size_t>,float>;

} // end of namespace


