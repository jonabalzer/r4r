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

#ifndef R4RTRAFO_H_
#define R4RTRAFO_H_

#include <vector>

#include "types.h"

namespace R4R {

/*! \brief transformation
 *
 * \todo Introduce template parameter for data type.
 *
 */
template<typename T,u_int n>
class CTransformation {

public:

	//! Constructor.
	CTransformation();

    //! Concatenates two transformations.
    CTransformation<T,n> operator*(const CTransformation<T,n>& x) const;

    //! Constructor.
    CTransformation(const CDenseArray<T>& F);

    //! Non-destructive access.
    T Get(u_int i, u_int j) const;

    //! Destructive access.
    T& operator()(u_int i, u_int j);

    //! Forward transformation.
    CVector<T,n> Transform(const CVector<T,n>& x);

    //! Forward transformation.
    CVector<T,n> Transform(const T* x);

    //! Typecast.
    operator CDenseArray<T>() const;

    //! Inverts the transformation.
    virtual bool Invert();

    //! Writes transformation to a stream.
    template <class U,u_int m> friend std::ostream& operator << (std::ostream& os, const CTransformation<U,m>& x);

    //! Reads transformation matrix from a stream.
    template <class U,u_int m> friend std::istream& operator >> (std::istream& is, CTransformation<U,m>& x);

    //! Low-level access to data.
    T* Data() { return m_F; }

    //! Parallelized mass transformation.
    std::vector<CVector<T,n> > Transform(const std::vector<CVector<T,n> >& x);

    //! Returns Jacobian of the transformation, i.e., its linear part.
    CDenseArray<T> GetJacobian();

    //! Access to translation vector.
    CVector<T,n> GetTranslation();

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

/*! \brief interface for infinitesimal rotations
 *
 *
 *
 */
template<typename T,u_int n>
class CDifferentialRotation:public CTransformation<T,n> {

public:

    //! Standard constructor.
    CDifferentialRotation():CTransformation<T,n>() {}

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

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

protected:

    using CTransformation<T,2>::m_F;

};

/*! \brief differential rotation in 2D
 *
 *
 *
 */
template<typename T>
class CDifferentialRotation<T,2>:public CTransformation<T,2> {

public:

    //! Constructor.
    CDifferentialRotation(T o);

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

    //! Logarithm.
    static void Log(const T* R, T& o1, T& o2, T& o3);

protected:

    using CTransformation<T,3>::m_F;

};


/*! \brief infinitesimal 3D rotation
 *
 *
 *
 */
template<typename T>
class CDifferentialRotation<T,3>:public CTransformation<T,3> {

public:

    //! Constructor.
    CDifferentialRotation(T o1, T o2, T o3, u_int dim);

    //! Rodrigues formula including derivatives.
    static void Rodrigues(const T& o1, const T& o2, const T& o3, T* R, T* DRo1, T* DRo2, T* DRo3);

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

    //! Constructor.
    CRigidMotion(const CDenseArray<T>& x):CTransformation<T,n>(x) {}

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
    CRigidMotion():CTransformation<T,2>() {}

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
    CRigidMotion():CTransformation<T,3>() {}

    //! Constructor.
    CRigidMotion(T o1, T o2, T o3, T t1, T t2, T t3);

    //! Constructor. FIXME: This has to project to SE(3)!
    CRigidMotion(const CDenseArray<T>& x):CTransformation<T,3>(x) {}

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

protected:

    using CTransformation<T,3>::m_F;

};



}



#endif /* TRAFO_H_ */
