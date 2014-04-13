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

#include "darray.h"

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

    //! Constructor which concatenates two transformations.
    CTransformation<T,n> operator*(const CTransformation<T,n>& x) const;

    //! Constructor.
    CTransformation(const CDenseArray<T>& F);

    //! Constructor using pointer to raw data.
    CTransformation(T* data);

    //! Non-destructive access.
    T Get(u_int i, u_int j) const;

    //! Destructive access.
    T& operator()(u_int i, u_int j);

    //! Forward transformation.
    CVector<T,n> Transform(const CVector<T,n>& x) const;

    //! Transform a direction.
    CVector<T,n> DifferentialTransform(const CVector<T,n>& x) const;

    //! Forward transformation.
    CVector<T,n> Transform(const T* x) const;

    //! Typecast.
    operator CDenseArray<T>() const;

    //! Inverts the transformation.
    virtual bool Invert();

    //! Writes transformation to a stream.
    template <class U,u_int m> friend std::ostream& operator << (std::ostream& os, const CTransformation<U,m>& x);

    //! Reads transformation matrix from a stream.
    template <class U,u_int m> friend std::istream& operator >> (std::istream& is, CTransformation<U,m>& x);

    //! Low-level access to data.
    const T* Data() const { return m_F; }

    //! Parallelized mass transformation.
    std::vector<CVector<T,n> > Transform(const std::vector<CVector<T,n> >& x) const;

    //! Parallelized mass transformation.
    std::vector<CVector<T,n> > DifferentialTransform(const std::vector<CVector<T,n> >& x) const;

    //! Returns Jacobian of the transformation, i.e., its linear part.
    CDenseArray<T> GetJacobian() const;

    //! Access to translation vector.
    CVector<T,n> GetTranslation() const;

    //! Access to the principal axis (e.g., of a pinhole camera).
    CVector<T,n> GetPrincipalAxis() const;

    //! Checks if two transformations are equal.
    bool operator==(const CTransformation<T,n>& x) const;


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

    //! Checks if two transformations are equal.
    bool operator==(CRigidMotion<T,n>& x) { return CTransformation<T,n>::operator ==(x); }

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
    CRigidMotion(T t1, T t2, T o);

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
    CRigidMotion(T t1, T t2, T t3, T o1, T o2, T o3);

    //! Constructor.
    CRigidMotion(const CVector<T,6>& m);

    //! Constructor. FIXME: This has to project to SE(3)!
    CRigidMotion(const CDenseArray<T>& x):CTransformation<T,3>(x) {}

    //! \copydoc CTransformation<T,n>::Invert()
    virtual bool Invert();

protected:

    using CTransformation<T,3>::m_F;

};

/*! \brief (time) sequence of coordinate frames
 *
 * The dimension is not a template parameter, because the excessive specialization
 * required would amount to implementing the class for each dimension separately
 * in the first place.
 *
 */
template<template<class U, class Allocator = std::allocator<U> > class Container,typename T>
class C3dTrajectory {

public:

    //! Constructor.
    C3dTrajectory():m_data() {}

    //! Inserts a new frame into the container.
    void Update(T t1, T t2, T t3, T o1, T o2, T o3);

    //! Inserts a new frame into the container.
    void Update(const CVector<T,6>& x) { m_data.push_back(x); }

    //! Inserts a new frame into the container.
    void Update(const CRigidMotion<T,3>& F);

    //! Returns latest state.
    CRigidMotion<T,3> GetLatestState() const;

    //! Returns the first state
    CRigidMotion<T,3> GetInitialState() const;

    //! Provides read-only access to the data.
    const Container<CVector<T,6> >& GetData() const { return m_data; }

    /*! Writes the trajectory to file.
     *
     * \param[in] filename file name
     * \param[in] format flag indicating whether to store the motion in
     * exponential coordinates or as full frames
     */
    bool WriteToFile(const char* filename, bool format = false) const;

private:

    Container<CVector<T,6> > m_data;

};


}



#endif /* TRAFO_H_ */
