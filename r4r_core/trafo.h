/*
 * trafo.h
 *
 *  Created on: Oct 16, 2012
 *      Author: jbalzer
 */

#ifndef R4RTRAFO_H_
#define R4RTRAFO_H_

#include "types.h"

namespace R4R {



/*! \brief transformation
 *
 *
 *
 */
template<typename T,u_int n>
class CTransformation {

public:

	//! Constructor.
	CTransformation();

    //! Non-destructive access.
    T Get(u_int i, u_int j) const;

    //! Destructive access.
    T& operator()(u_int i, u_int j);

    //! Forward transformation.
    CVector<T,n> Transform(const CVector<T,n>& x);

    //! Typecast.
    operator CDenseArray<T>() const;

    //! In-place forward transformation.
    void Transform(T* x);

    //! Inverts the transformation.
    virtual bool Invert();

    //! Writes transformation to a stream.
    template <class U,u_int m> friend std::ostream& operator << (std::ostream& os, const CTransformation<U,m>& x);

    //! Low-level access to data.
    T* Data() { return m_F; }

protected:

    T m_F[n*(n+1)];               //!< container for data

};

/*! \brief rotation interface
 *
 *
 *
 */
template<typename T,u_int n>
class CRotation:public CTransformation<T,n> {

public:

    //! Standard constructor.
    CRotation():CTransformation<T,n>() {}

protected:

    using CTransformation<T,n>::m_F;

};


/*! \brief 2D rotation
 *
 *
 *
 */
template<typename T>
class CRotation<T,2>:public CTransformation<T,2> {

public:

    //! Constructor.
    CRotation(T o);

    //! Rodrigues formula.
    static void Rodrigues(const T& o, T* R);

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

protected:

    using CTransformation<T,2>::m_F;

};

/*! \brief 3D rotation
 *
 *
 *
 */
template<typename T>
class CRotation<T,3>:public CTransformation<T,3> {

public:

    //! Constructor.
    CRotation(T o1, T o2, T o3);

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

    //! Rodrigues formula.
    static void Rodrigues(const T& o1, const T& o2, const T& o3, T* R);

    //! Rodrigues formula including derivatives.
    static void Rodrigues(const T& o1, const T& o2, const T& o3, T* R, T* DRo1, T* DRo2, T* DRo3);

    //! Logarithm.
    static void Log(const T* R, T& o1, T& o2, T& o3);

protected:

    using CTransformation<T,3>::m_F;

};

/*! \brief rigid motion interface
 *
 *
 *
 */
template<typename T,u_int n>
class CRigidMotion:public CTransformation<T,n> {

public:

    //! Constructor.
    CRigidMotion():CTransformation<T,n>() {}

protected:

    using CTransformation<T,n>::m_F;

};

/*! \brief rigid motion in 2D
 *
 *
 *
 */
template<typename T>
class CRigidMotion<T,2>:public CTransformation<T,2> {

public:

    //! Constructor.
    CRigidMotion(T o, T t1, T t2);

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

protected:

    using CTransformation<T,2>::m_F;

};

/*! \brief rigid motion in 3D
 *
 *
 *
 */
template<typename T>
class CRigidMotion<T,3>:public CTransformation<T,3> {

public:

    //! Constructor.
    CRigidMotion(T o1, T o2, T o3, T t1, T t2, T t3);

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

protected:

    using CTransformation<T,3>::m_F;

};



}



#endif /* TRAFO_H_ */
