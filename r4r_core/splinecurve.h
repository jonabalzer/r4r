#ifndef CSPLINECURVE_H
#define CSPLINECURVE_H


#include "darray.h"

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
    CSplineCurve(u_int d, u_int p, u_int n);

    //! Constructor.
    CSplineCurve(const CDenseVector<T>& knot, const CDenseArray<T>& cv);

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

    u_int m_d;                              //!< dimension of the embedding space
    u_int m_p;                              //!< degree
    u_int m_n;                              //!< number of control points
    u_int m_k;                              //!< number of knots
    CDenseVector<T> m_knot;                 //!< storage for knot vector
    CDenseArray<T> m_cv;                    //!< storage for control points


};

} // end namespace

#endif // CSPLINECURVE_H
