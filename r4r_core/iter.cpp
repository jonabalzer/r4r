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

#include "iter.h"
#include "darray.h"
#include "rutils.h"

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>

namespace R4R {

using namespace std;

template<class Matrix,typename T>
CConjugateGradientMethod<Matrix,T>::CConjugateGradientMethod(const CPreconditioner<Matrix,T>& M, size_t n, double eps, bool silent):
    CIterativeLinearSolver<Matrix,T>::CIterativeLinearSolver(M,n,eps,silent) {}

template<class Matrix,typename T>
vector<double> CConjugateGradientMethod<Matrix,T>::Iterate(const Matrix& A, const CDenseArray<T>& B, CDenseArray<T>& X) const {

    // check dimensions
    if(!(A.NCols()==X.NRows() && X.NRows()==B.NRows() && X.NCols()==B.NCols())) {

        cerr << "ERROR: Check matrix dimensions!" << endl;
        return vector<double>();

    }

    // init iteration index
    size_t k = 0;

    // compute residual
    CDenseArray<T> R = B - A*X;

    // residual norm
    vector<double> res;
    res.push_back(R.Norm2());

    if(!m_silent)
        cout << "k=" << k << ": " << res.back() << endl;

    // make a copy for the temporary variable used in pre-conditioning
    CDenseArray<T> Z = R.Clone();

    // apply preconditioner
    m_M.Solve(Z,R);

    CDenseVector<T> deltan = CDenseArray<T>::ColumwiseInnerProduct(Z,R);

    // init descent direction
    CDenseArray<T> P = Z.Clone();

    while(k<m_n) {

        CDenseArray<T> Q = A*P;

        CDenseVector<T> pq =  CDenseArray<T>::ColumwiseInnerProduct(P,Q);
        CDenseVector<T> alpha = deltan/pq;

        // descent step
        X = X + P.ScaleColumns(alpha);

        // update residual
        Q = Q.ScaleColumns(alpha);
        R = R - Q;

        // check convergence
        res.push_back(R.Norm2());

        k++;

        if(!m_silent)
            cout << "k=" << k << ": " << res.back() << endl;

        if(res.back()<m_eps)
            break;

        // apply pre-conditioner
        m_M.Solve(Z,R);

        // make a copy of deltan
        CDenseVector<T> deltao = deltan.Clone();

        // column-wise inner product
        deltan = CDenseVector<T>::ColumwiseInnerProduct(Z,R);

        // component-wise division
        CDenseVector<T> beta = deltan/deltao;

        // update descent direction
        P = Z + P.ScaleColumns(beta);

    }

    return res;

}

template<class Matrix,typename T>
vector<double> CConjugateGradientMethod<Matrix,T>::Iterate(const Matrix& A, const CDenseVector<T>& b, CDenseVector<T>& x) const {

    // check dimensions
    if(!(A.NCols()==x.NRows() && x.NRows()==b.NRows())) {

        cerr << "ERROR: Check matrix dimensions!" << endl;
        return vector<double>();

    }

    // init
    size_t k = 0;

    CDenseVector<T> r = b - A*x;

    vector<double> res;
    res.push_back(r.Norm2());

    if(!m_silent)
        cout << "k=" << k << ": " << res.back() << endl;

    CDenseVector<T> z = r.Clone();

    m_M.Solve(z,r);

    T deltao = CDenseArray<T>::InnerProduct(z,r);

    CDenseVector<T> p = z.Clone();

    while(k<m_n) {

        CDenseVector<T> q = A*p;

        T alpha = deltao/CDenseArray<T>::InnerProduct(p,q);

        x = x + p*alpha;

        q.Scale(alpha);
        r = r - q;

        res.push_back(r.Norm2());

        k++;

        if(!m_silent)
            cout << "k=" << k << ": " << res.back() << endl;

        if(res.back()<m_eps)
            break;


        m_M.Solve(z,r);

        T deltan = CDenseArray<T>::InnerProduct(z,r);

        T beta = deltan/deltao;

        p.Scale(beta);

        p = z + p;

        deltao = deltan;

    }

    return res;

}

template class CConjugateGradientMethod<CDenseArray<double>,double>;
template class CConjugateGradientMethod<CSparseArray<double>,double>;
template class CConjugateGradientMethod<CDenseArray<float>,float>;
template class CConjugateGradientMethod<CSparseArray<float>,float>;

template<class Matrix,typename T>
CConjugateGradientMethodLeastSquares<Matrix,T>::CConjugateGradientMethodLeastSquares(const CPreconditioner<Matrix,T>& M, size_t n, double eps, bool silent):
    CIterativeLinearSolver<Matrix,T>::CIterativeLinearSolver(M,n,eps,silent) {}

template<class Matrix,typename T>
vector<double> CConjugateGradientMethodLeastSquares<Matrix,T>::Iterate(const Matrix& A, const CDenseArray<T>& B, CDenseArray<T>& X) const {

    if(!(A.NCols()==X.NRows() && A.NRows()==B.NRows() && X.NCols()==B.NCols())) {

        cerr << "ERROR: Check matrix dimensions!" << endl;
        return vector<double>();

    }

    /* Get a transposed copy of the input matrix. This way we can keep the input reference
     * constant (would not work for in-place back and forth transposition). The overhead is
     * minimal for the dense matrix structure because smart pointers are used. Need to introduce
     * same structure to sparse matrix classes.
     *
     */
    Matrix At = Matrix::Transpose(A);

    // init
    size_t k = 0;

    // residual of the non-square system
    CDenseArray<T> R = B - A*X;

    // strore it
    vector<double> res;
    res.push_back(R.Norm2());

    // residual of the normal equation
    CDenseArray<T> Rnormal = At*R;

    // preconditioning
    CDenseArray<T> Z = Rnormal.Clone();
    m_M.Solve(Z,Rnormal);

    CDenseVector<T> deltao = CDenseArray<T>::ColumwiseInnerProduct(Z,Rnormal);

    // descent direction
    CDenseArray<T> P = Z.Clone();

    if(!m_silent)
        cout << "k=" << k << ": " << res.back() << endl;

    if(res.back()<m_eps)
        return res;

    while(k<m_n) {

        // need that later
        CDenseArray<T> Q = A*P;

        // step siz
        CDenseVector<T> qq = CDenseArray<T>::ColumwiseInnerProduct(Q,Q);
        CDenseVector<T> alpha = deltao/qq;

        // perform descent step
        X = X + P.ScaleColumns(alpha);

        // update residual of non-square system
        R = R - Q.ScaleColumns(alpha);

        // check convergence
        res.push_back(R.Norm2());

        k++;

        if(!m_silent)
            cout << "k=" << k << ": " << res.back() << endl;

        // FIXME: add another constant that monitors absolute value of residual
        if(fabs(res.at(res.size()-2)-res.back())<m_eps || res.back()<m_eps)
            break;

        // update residual of normal equation
        Rnormal = At*R;

        // apply preconditioner
        m_M.Solve(Z,Rnormal);

        // update beta
        CDenseVector<T> deltan = CDenseArray<T>::ColumwiseInnerProduct(Z,Rnormal);
        CDenseVector<T> beta = deltan/deltao;
        deltao = deltan;

        // update direction
        P = Z + P.ScaleColumns(beta);

    }

    return res;

}

template<class Matrix,typename T>
vector<double> CConjugateGradientMethodLeastSquares<Matrix,T>::Iterate(const Matrix& A, const CDenseVector<T>& b, CDenseVector<T>& x) const {

    if(!(A.NCols()==x.NRows() && A.NRows()==b.NRows())) {

        cerr << "ERROR: Check matrix dimensions!" << endl;
        return vector<double>();

    }

    /* Get a transposed copy of the input matrix. This way we can keep the input reference
     * constant (would not work for in-place back and forth transposition). The overhead is
     * minimal for the dense matrix structure because smart pointers are used. Need to introduce
     * same structure to sparse matrix classes.
     *
     */
    Matrix At = Matrix::Transpose(A);

    // init
    size_t k = 0;

    // residual of the non-square system
    CDenseVector<T> r = b - A*x;

    // residuals
    vector<double> res;
    res.push_back(r.Norm2());

    // residual of the normal equation
    //A.Transpose();
    CDenseVector<T> rnormal = At*r;
    //A.Transpose();

    // preconditioning
    CDenseVector<T> z = rnormal.Clone();
    cout << "FICK MICH" << z << endl;
    cout << "FICK MICH AUCH" << rnormal << endl;

    m_M.Solve(z,rnormal);


    // descent direction
    CDenseVector<T> p = z.Clone();

    T deltao = CDenseArray<T>::InnerProduct(z,rnormal);

    if(!m_silent)
        cout << "k=" << k << ": " << res.back() << endl;

    while(k<m_n) {

        // need that later
        CDenseVector<T> q = A*p;

        // step size (recycle deltao for computation of beta)
        T alpha = deltao/CDenseArray<T>::InnerProduct(q,q);

        // perform descent step
        x = x + p*alpha;

        q.Scale(alpha);
        r = r - q;				// update residual of non-square system

        res.push_back(r.Norm2());

        k++;

        if(!m_silent)
            cout << "k=" << k << ": " << res.back() << endl;

        if(fabs(res.at(res.size()-2)-res.back())<m_eps)
            break;

        // update residual of normal equation
        rnormal = At*r;

        // apply preconditioner
        m_M.Solve(z,rnormal);

        // update beta
        T deltan = CDenseArray<T>::InnerProduct(z,rnormal);
        T beta = deltan/deltao;
        deltao = deltan;

        // update direction
        p.Scale(beta);
        p = z + p;

    }

    return res;

}

template class CConjugateGradientMethodLeastSquares<CDenseArray<double>,double>;
template class CConjugateGradientMethodLeastSquares<CSparseArray<double>,double>;
template class CConjugateGradientMethodLeastSquares<CDenseArray<float>,float>;
template class CConjugateGradientMethodLeastSquares<CSparseArray<float>,float>;

}
