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

/*! \brief iterative linear solver interface
 *
 *
 *
 */
template<class Matrix,typename T>
class CIterativeLinearSolver {

protected:

    /*! \brief Constructor.
     * \param[in] M preconditioner
     * \param[in] n (maximum) number of iterations
     * \param[in] eps absolute residual at which to terminate
     * \param[in] silent verbosity flag
     *
     * This class has pure virtual functions. Block direct creation.
     *
     */
    CIterativeLinearSolver(const CPreconditioner<Matrix,T>& M, size_t n, double eps, bool silent = true):m_M(M),m_n(n),m_eps(eps),m_silent(silent) {}

    //! Standard constructor (deleted).
    CIterativeLinearSolver() = delete;

public:

    //! Access to accuracy.
    void SetEps(double eps) { m_eps = eps; }

    //! Access to maximum number of iterations.
    void SetN(size_t n) { m_n = n; }

    //! Set silent.
    void SetSilent(bool silent) { m_silent = silent; }

    /*! \brief Iterate.
     *
     * \param[in] A matrix \f$A\in\mathbb{R}^{m\times n}\f$
     * \param[in] B right-hand side \f$B\in\mathbb{R}^{m\times d}\f$
     * \param[in] X solution \f$X\in\mathbb{R}^{n\times d}\f$
     * \returns residual norm
     *
     * This function solves the system for \f$d\f$ right-hand sides organized as columns
     * of a matrix \f$B\in\mathbb{R}^{m\times d}\f$. The iteration is started from the
     * value of \f$X\f$ which is passed to the function.
     *
     */
    virtual double Iterate(const Matrix& A, const CDenseArray<T>& b, CDenseArray<T>& x) = 0;

    /*! \brief Iterate.
     *
     * \param[in] A matrix \f$A\in\mathbb{R}^{m\times n}\f$
     * \param[in] b right-hand side \f$b\in\mathbb{R}^{m}\f$
     * \param[in] x solution \f$x\in\mathbb{R}^{n}\f$
     * \returns residual norm
     *
     * This function solves the system for a single right-hand side. The implementation is
     * optimized for this special case.
     *
     */
    virtual double Iterate(Matrix& A, const CDenseVector<T>& b, CDenseVector<T>& x) = 0;

protected:

    const CPreconditioner<Matrix,T>& m_M;    		//!< preconditioner
    size_t m_n;                             		//!< number of steps
    double m_eps;                               	//!< absolute accuracy
    bool m_silent;                                  //!< flag that determines whether messages are displayed

};

/*! \brief Conjugate gradient method.
 *
 * Implements the conjugate gradient (CG) method, cf. [Hestenes1952]. It is crucial that the matrix \f$A\f$
 * be symmetric and positive-definite.
 *
 */
template<class Matrix,typename T>
class CConjugateGradientMethod: public CIterativeLinearSolver<Matrix,T> {

public:

    //! \copybrief CIterativeLinearSolver::CIterativeLinearSolver(CPreconditioner&,size_t,T,bool)
    CConjugateGradientMethod(const CPreconditioner<Matrix,T>& M, size_t n, double eps, bool silent = true);
    CConjugateGradientMethod() = delete;

    //! \copydoc CIterativeLinearSolver::Iterate(const Matrix&,const CDenseArray<T>&,CDenseArray<T>&)
    double Iterate(const Matrix& A, const CDenseArray<T>& B, CDenseArray<T>& X);

    //! \copydoc CIterativeLinearSolver::Iterate(const Matrix&,const CDenseVector<T>&,CDenseVector<T>&)
    double Iterate(Matrix& A, const CDenseVector<T>& b, CDenseVector<T>& x);

private:

    using CIterativeLinearSolver<Matrix,T>::m_M;
    using CIterativeLinearSolver<Matrix,T>::m_n;
    using CIterativeLinearSolver<Matrix,T>::m_eps;
    using CIterativeLinearSolver<Matrix,T>::m_silent;

};


/*! \brief Conjugate gradient least-squares method.
 *
 * Implements the conjugate gradient least-squares (CGLS) method, cf. [Paige1982]. The linear system
 * may be over-determined. A solution is obtained without explicit formation of the normal equation.
 *
 */
template<class Matrix,typename T>
class CConjugateGradientMethodLeastSquares:public CIterativeLinearSolver<Matrix,T> {

public:

    //! \copybrief CIterativeLinearSolver::CIterativeLinearSolver(CPreconditioner&,size_t,T,bool)
    CConjugateGradientMethodLeastSquares(const CPreconditioner<Matrix,T>& M, size_t n, double eps, bool silent = true);
    CConjugateGradientMethodLeastSquares() = delete;

    //! \copydoc CIterativeLinearSolver::Iterate(const Matrix&,const CDenseArray<T>&,CDenseArray<T>&)
    double Iterate(const Matrix& A, const CDenseArray<T>& B, CDenseArray<T>& X);

    //! \copydoc CIterativeLinearSolver::Iterate(const Matrix&,const CDenseVector<T>&,CDenseVector<T>&)
    double Iterate(Matrix& A, const CDenseVector<T>& b, CDenseVector<T>& x);

private:

    using CIterativeLinearSolver<Matrix,T>::m_M;
    using CIterativeLinearSolver<Matrix,T>::m_n;
    using CIterativeLinearSolver<Matrix,T>::m_eps;
    using CIterativeLinearSolver<Matrix,T>::m_silent;

};


}

#endif /* ITER_H_ */
