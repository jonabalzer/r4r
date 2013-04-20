#ifndef CSPLINECURVE_H
#define CSPLINECURVE_H


#include <memory>
#include <xmmintrin.h>
#include <malloc.h>

class  CSplineDeallocator {

public:

    void operator()(float* p) const { _mm_free(p); }

};

class CSplineCurve
{

public:

    //! Constructor.
    CSplineCurve();

    //! Constructor.
    CSplineCurve(size_t d, size_t p, size_t n);

    //! Destructor.
    virtual ~CSplineCurve();

    /*! \brief Makes uniform clamped knot vector.
     *
     * Clamped means there a \f$p\f$ knots at the end points (due to
     * strange ON convention.
     *
     */
    void MakeClampedUniformKnotVector(float a, float b);

    /*! \brief Returns the knot span.
     *
     * \details How should we do this? If the knot vector is uniform, this can
     * be done efficiently. If not, we need to run a binary search.
     *
     */
    int GetSpan(float t);

    //! Access to control points.
    std::shared_ptr<float> GetCVData() { return m_cv; }

    //! Access to knot vector.
    std::shared_ptr<float> GetKnotVector() { return m_knot; }

    //! Prints information about the spline.
    void Print();

    //! Cox-de-Boor formula.
    static bool EvaluateNurbsBasis(int order, const float* knot, float t, float* N);

    //! Cox-de-Boor recursion computing derivatives.
    static bool EvaluateNurbsBasisDerivatives(int order, const float* knot, int der_count, float* N);

private:

    size_t m_d;                         //!< dimension of the embedding space
    size_t m_p;                         //!< degree
    size_t m_n;                         //!< number of control points
    size_t m_k;                         //!< number of knots
    std::shared_ptr<float> m_knot;           //!< storage for knot vector
    std::shared_ptr<float> m_cv;             //!< storage for control points

    // protect assignment and copy
    CSplineCurve(const CSplineCurve& x);
    CSplineCurve& operator =(const CSplineCurve& x);

};

#endif // CSPLINECURVE_H
