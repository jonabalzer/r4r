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

#ifndef R4RVECN_H
#define R4RVECN_H

#include <stdlib.h>
#include <iostream>

namespace R4R {

template<typename T> class CDenseVector;

/*! \brief short vectors of length \f$n\f$
 *
 *
 *
 */
template<typename T,u_int n>
class CVector {

public:

    //! Constructor.
    CVector();

    //! Constructor.
    CVector(T val);

    //! Constructor.
    CVector(const CDenseVector<T>& x);

    //! Initializer list constructor.
    CVector(std::initializer_list<T> list);

    //! Typecast operator.
    template<typename U> operator CVector<U,n>() const {

        CVector<U,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = U(Get(i));

        return result;

    }

    //! Computes \f$l_p\f$-norm.
    double Norm(double p) const;

    //! Computes \f$l_2\f$-norm.
    double Norm2() const;

    //! Normalizes the vector.
    bool Normalize();

    //! Read-write element access.
    T& operator()(u_int i);

    //! In-place addition.
    void operator+=(const CVector<T,n>& x);

    //! In-place element-wise multiplication.
    void operator*=(const CVector<T,n>& x);

    //! Read element access.
    T Get(u_int i) const;

    //! Writes vector to a stream.
    template <class U,u_int m> friend std::ostream& operator << (std::ostream& os, const CVector<U,m>& x);

    //! Reads vector from a stream.
    template <class U,u_int m> friend std::istream& operator >> (std::istream& is, CVector<U,m>& x);

    //! Low-level acces to the data.
    T* Data() { return m_data; }

private:

    T m_data[n];

};

template <typename T>
inline CVector<T,3> Cross(const CVector<T,3>& x, const CVector<T,3>& y) {

    CVector<T,3> result;

    result(0) = x.Get(1)*y.Get(2) - x.Get(2)*y.Get(1);
    result(1) = x.Get(2)*y.Get(0) - x.Get(0)*y.Get(2);
    result(2) = x.Get(0)*y.Get(1) - x.Get(1)*y.Get(0);

    return result;

}

//! Dot product between to vectors.
template <typename T,u_int n>
inline double InnerProduct(const CVector<T,n>& x, const CVector<T,n>& y) {

    double result = 0;

    for(u_int i=0; i<n; i++)
        result += (double)x.Get(i)*(double)y.Get(i);

    return result;

}


//! Adds two vectors.
template <typename T,u_int n>
inline CVector<T,n> operator+(const CVector<T,n>& x, const CVector<T,n>& y) {

      CVector<T,n> result;

      for(u_int i=0; i<n; i++)
          result(i) = x.Get(i) + y.Get(i);

      return result;

}

//! Subtracts two vectors.
template <typename T,u_int n>
inline CVector<T,n> operator-(const CVector<T,n>& x, const CVector<T,n>& y) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i) - y.Get(i);

    return result;

}

//! Multiplies two vectors element-wise.
template <typename T,u_int n>
inline CVector<T,n> operator*(const CVector<T,n>& x, const CVector<T,n>& y) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i)*y.Get(i);

    return result;

}

//! Divides two vectors element-wise.
template <typename T,u_int n>
inline CVector<T,n> operator/(const CVector<T,n>& x, const CVector<T,n>& y) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i)/y.Get(i);

    return result;

}

//! Post-multiplies a vector by a scalar.
template <typename T,u_int n, typename U>
inline CVector<T,n> operator*(const CVector<T,n>& x, const U& s) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i)*(T)s;

    return result;

}

//! Pre-multiplies a vector by a scalar.
template <typename T,u_int n,typename U>
CVector<T,n> operator*(const U& s, const CVector<T,n>& x);

//! Post-divides a vector by a scalar.
template <typename T,u_int n, typename U>
inline CVector<T,n> operator/(const CVector<T,n>& x, const U& s) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i)/(T)s;

    return result;

}

//! Post-adds a scalar to the vector.
template <typename T,u_int n, typename U>
inline CVector<T,n> operator+(const CVector<T,n>& x, const U& s) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i) + (T)s;

    return result;

}

//! Pre-adds a scalar to the vector.
template <typename T,u_int n,typename U>
inline CVector<T,n> operator+(const U& s, const CVector<T,n>& x) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i) + (T)s;

    return result;

}

//! Post-subtracts a scalar from the vector.
template <typename T,u_int n, typename U>
inline CVector<T,n> operator-(const CVector<T,n>& x, const U& s) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = x.Get(i) - (T)s;

    return result;

}

//! Pre-subtracts a scalar from the vector.
template <typename T,u_int n,typename U>
inline CVector<T,n> operator-(const U& s, const CVector<T,n>& x) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = (T)s - x.Get(i);

    return result;

}

//! Checks two vectors for equality.
template <typename T,u_int n>
inline bool operator==(const CVector<T,n>& x, const CVector<T,n>& y) {

    bool result = true;

    for(u_int i=0; i<n; i++)
        result = result && (x.Get(i)==y.Get(i));

    return result;

}

//! Checks two vectors for inequality.
template <typename T,u_int n>
inline bool operator!=(const CVector<T,n>& x, const CVector<T,n>& y) {

    return !(x==y);

}


}



#endif // VECN_H
