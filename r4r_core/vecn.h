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

#ifndef R4RVECN_H
#define R4RVECN_H

#include <stdlib.h>
#include <iostream>

#ifdef HAVE_EXR
#include <half.h>
#endif


namespace R4R {

template<typename T> class CDenseVector;

enum class DIM {  ZERO = 0, ONE = 1, TWO = 2, THREE = 3 };

/*! \brief short vectors of length \f$n\f$
 *
 *
 *
 */
template<typename T,u_int n>
class CVector {

private:

    T m_data[n];

public:

    //! Constructor.
    CVector();

    //! Constructor.
    CVector(T val);

    //! Constructor.
    CVector(const CDenseVector<T>& x);

    //! Initializer list constructor.
    CVector(std::initializer_list<T> list);

    //! Fast routine for checking whether the vector is zero.
    bool IsZero() const;

    //! Computes \f$l_p\f$-norm.
    double Norm(double p) const;

    //! Computes \f$l_2\f$-norm.
    double Norm2() const;

    //! Computes the absolute value.
    CVector<T,n> Abs() const;

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
    const T* Data() const { return m_data; }

    //! Maximum.
    T Max();

    //! Typecast operator.
    template<typename U> operator CVector<U,n>() const {

        CVector<U,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = U(Get(i));

        return result;

    }

    //! Elementwise maximum.
    friend CVector<T,n> Max(const CVector<T,n>& x, const CVector<T,n>& y) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = std::max<T>(x.Get(i),y.Get(i));

        return result;

    }

    //! Elementwise minimum.
    friend CVector<T,n> Min(const CVector<T,n>& x, const CVector<T,n>& y) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = std::min<T>(x.Get(i),y.Get(i));

        return result;

    }

    //! Adds two vectors.
    friend CVector<T,n> operator+(const CVector<T,n>& x, const CVector<T,n>& y) {

          CVector<T,n> result;

          for(u_int i=0; i<n; i++)
              result(i) = x.Get(i) + y.Get(i);

          return result;
    }

    /*! Dot product between to vectors.
     *
     * \TODO Add template parameter for return type.
     *
     */
    friend T InnerProduct(const CVector<T,n>& x, const CVector<T,n>& y) {

        T result = 0;

        for(u_int i=0; i<n; i++)
            result += x.Get(i)*y.Get(i);

        return result;

    }

    //! Cross-product between to vectors
    friend CVector<T,n> Cross(const CVector<T,n>& x, const CVector<T,n>& y) {

        // only defined for n=3
        if(n==3) {

            CVector<T,3> result;

            result(0) = x.Get(1)*y.Get(2) - x.Get(2)*y.Get(1);
            result(1) = x.Get(2)*y.Get(0) - x.Get(0)*y.Get(2);
            result(2) = x.Get(0)*y.Get(1) - x.Get(1)*y.Get(0);

            return result;

        }

        return CVector<T,n>();

    }

    //! Subtracts two vectors.
    friend CVector<T,n> operator-(const CVector<T,n>& x, const CVector<T,n>& y) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i) - y.Get(i);

        return result;

    }

    //! Multiplies two vectors element-wise.
    friend CVector<T,n> operator*(const CVector<T,n>& x, const CVector<T,n>& y) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i)*y.Get(i);

        return result;

    }

    //! Divides two vectors element-wise.
    friend CVector<T,n> operator/(const CVector<T,n>& x, const CVector<T,n>& y) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i)/y.Get(i);

        return result;

    }

    /*! Post-multiplies a vector by a scalar.
     *
     */
    template <typename U>
    friend CVector<T,n> operator*(const CVector<T,n>& x, const U& s) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i)*(T)s;

        return result;

    }

    //! Post-divides a vector by a scalar.
    template <typename U>
    friend CVector<T,n> operator/(const CVector<T,n>& x, const U& s) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i)/(T)s;

        return result;

    }


    //! Post-adds a scalar to the vector.
    template <typename U>
    friend CVector<T,n> operator+(const CVector<T,n>& x, const U& s) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i) + (T)s;

        return result;

    }

    //! Pre-adds a scalar to the vector.
    template <typename U>
    friend CVector<T,n> operator+(const U& s, const CVector<T,n>& x) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i) + (T)s;

        return result;

    }

    //! Post-subtracts a scalar from the vector.
    template<typename U>
    friend CVector<T,n> operator-(const CVector<T,n>& x, const U& s) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = x.Get(i) - (T)s;

        return result;

    }

    //! Pre-subtracts a scalar from the vector.
    template<typename U>
    friend CVector<T,n> operator-(const U& s, const CVector<T,n>& x) {

        CVector<T,n> result;

        for(u_int i=0; i<n; i++)
            result(i) = (T)s - x.Get(i);

        return result;

    }

    //! Checks two vectors for equality.
    friend bool operator==(const CVector<T,n>& x, const CVector<T,n>& y) {

        bool result = true;

        for(u_int i=0; i<n; i++)
            result = result && (x.Get(i)==y.Get(i));

        return result;

    }

    //! Checks two vectors for inequality.
    friend bool operator!=(const CVector<T,n>& x, const CVector<T,n>& y) {

        return !(x==y);

    }

    /*! Checks two vectors for inequality.
     *
     * This is default order relation is induced by the positive cone. Use comparator
     * class for a lexicographical order relation.
     *
     */
    friend bool operator<(const CVector<T,n>& x, const CVector<T,n>& y) {

        bool result = true;

        for(u_int i=0; i<n; i++)
            result = result && (x.Get(i)<y.Get(i));

        return result;

    }

    /*! Checks two vectors for inequality.
     *
     * This is default order relation is induced by the positive cone. Use comparator
     * class for a lexicographical order relation.
     *
     */
    friend bool operator<=(const CVector<T,n>& x, const CVector<T,n>& y) {

        bool result = true;

        for(u_int i=0; i<n; i++)
            result = result && (x.Get(i)<=y.Get(i));

        return result;

    }

};


template<typename T,u_int n>
class CConeOrder {

    //! Compares two views based on their index.
    bool operator()(const CVector<T,n>& x, const CVector<T,n>& y) {

        bool result = true;

        for(u_int i=0; i<n; i++)
            result = result && (x.Get(i)<=y.Get(i));

        return result;

    }

};

/*! Pre-multiplies a vector by a scalar.
 *
 * This is instantiated explicitly to avoid conflicts with matrix-vector
 * multiplication of #CDenseArray.
 *
 */
template <typename T,u_int n,typename U>
CVector<T,n> operator*(const U& s, const CVector<T,n>& x);

template<typename T,u_int n>
class CLexicographicOrder {

public:

    //! Lexicographic order relation.
    bool operator()(const CVector<T,n>& x, const CVector<T,n>& y) {

        // start recursion
        return this->Compare(x,y,0);

    }

private:

    bool Compare(const CVector<T,n>& x, const CVector<T,n>& y, u_int i) {

        // components are equal?
        if(i==n-1)
            return x.Get(i)<=y.Get(i);

        if (x.Get(i)>y.Get(i))
            return false;
        else
            return Compare(x,y,i+1);

    }

};


typedef CVector<double,3> vec3;
typedef CVector<float,3> vec3f;
typedef CVector<double,2> vec2;
typedef CVector<float,2> vec2f;
typedef CVector<unsigned char,3> rgb;

#ifdef HAVE_EXR
typedef CVector<half,3> vec3h;
#endif

}



#endif // VECN_H
