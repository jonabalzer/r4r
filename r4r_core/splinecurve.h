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

#ifndef R4RSPLINECURVE_H
#define R4RSPLINECURVE_H


#include "darray.h"
#include "kernels.h"

namespace R4R {

/*! \brief B-spline curve in \f$\mathbb{R}^d\f$
*
*
*
*
*
*
*/
template <class T>
class CSplineCurve
{

public:

    //! Constructor.
    CSplineCurve();

    //! Constructor.
    CSplineCurve(size_t d, size_t p, size_t n);

    //! Constructor.
    CSplineCurve(const CDenseVector<T>& knot, const CDenseArray<T>& cv);

    /*! \brief Constructor.
     *
     * \details This constructor reconstructs the spline parameters from the set of control
     * points. It assumes a uniform clamped knot vector supported in \f$[0,1]\f$.
     *
     */
    CSplineCurve(const CDenseArray<T>& cv, size_t p = 2);

    /*! \brief Makes uniform clamped knot vector.
     *
     * Clamped means there a \f$p\f$ knots at the end points (due to
     * strange ON convention.
     *
     */
    void MakeClampedUniformKnotVector(T a, T b);

    /*! \brief Returns the knot span.
     *
     * \details How should we do this? If the knot vector is uniform, this can
     * be done efficiently. If not, we need to run a binary search.
     *
     */
    int GetSpan(T t);

    //! Finds the closest point on the curve to a given point \f$x\f$.
    CDenseVector<T> FindLocallyClosestPoint(const CDenseVector<T>& y, CMercerKernel<T>& kernel, const T hint, const T eps);

    //! Evaluates the spline.
    CDenseVector<T> Evaluate(T t);

    //! Computes the tangent to the spline curve.
    CDenseArray<T> Tangent(T t);

    //! Computes the curvature vector.
    CDenseArray<T> Normal(T t);

    //! Distance of a point to the curve.
    T Distance(const CDenseVector<T>& x);

    //! Distance of another spline curve to the curve.
    T Distance(const CSplineCurve<T>& curve);

    //! Applies a rigid transformation to the control points.
    void Transform(const CDenseArray<T>& M);

    //! Access to control points.
    CDenseArray<T>& GetCVData() { return m_cv; }

    //! Access to knot vector.
    CDenseVector<T>& GetKnotVector() { return m_knot; }

    //! Prints information about the spline.
    void Print();

    //! Cox-de-Boor formula.
    static bool EvaluateNurbsBasis(u_int order, const T* knot, T t, T* N);

    //! Cox-de-Boor recursion computing derivatives.
    static bool EvaluateNurbsBasisDerivatives(u_int order, const T* knot, u_int der_count, T* N);

private:

    size_t m_d;                              //!< dimension of the embedding space
    size_t m_p;                              //!< degree
    size_t m_n;                              //!< number of control points
    size_t m_k;                              //!< number of knots
    CDenseVector<T> m_knot;                  //!< storage for knot vector
    CDenseArray<T> m_cv;                     //!< storage for control points


};

} // end namespace

#endif // CSPLINECURVE_H
