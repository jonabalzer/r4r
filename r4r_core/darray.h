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

#ifndef R4RDARRAY_H_
#define R4RDARRAY_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <memory>

#ifdef __SSE4_1__
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <malloc.h>
#endif

#include "kernels.h"
#include "vecn.h"

namespace R4R {

template<class T>
class  CDenseMatrixDeallocator {

public:

#ifndef __SSE4_1__
    void operator()(T* p) const { delete [] p; }
#else
    void operator()(T* p) const { _mm_free(p); }
#endif

};

// use typesafe enums in C++11
enum class ETYPE {  NA = 0,
                    B1U = 1,
                    C1U = 2, C1S = 3,
                    S2U = 4, S2S = 5,
                    I4S = 6, I4U = 7,
                    F4S = 8,
                    L8S = 9, L8U = 10,
                    D8S = 11,
                    C3U = 12 };

template<class T> class CDenseVector;

/*! \brief dense 2d matrix/array
 *
 *
 *
 */
template<class T>
class CDenseArray {

public:

	//! Constructor
	CDenseArray();

	//! Constructor
    CDenseArray(size_t nrows, size_t ncols, T val = T(0));

    //! Re-allocation.
    void Resize(size_t nrows, size_t ncols);

	//! Copy constructor.
	CDenseArray(const CDenseArray& array);

    //! Concatenate two arrays.
    void Concatenate(const CDenseArray& array, bool direction);

    //! Assignment operator.
    CDenseArray<T> operator=(const CDenseArray<T>& array);

    //! Copy constructor.
    CDenseArray(size_t nrows, size_t ncols, std::shared_ptr<T> data);

    //! Hardcopy.
    CDenseArray<T> Clone();

	//! Destructor.
	virtual ~CDenseArray();

	//! Creates an identity matrix.
	void Eye();

	//! Fills in ones.
	void Ones();

	//! Files the array with uniformly distributed random numbers between \f$0\f$ and \f$1\f$.
    void Rand(T min, T max);

    //! Files the array with normally distributed random numbers.
    void RandN(T mu, T sigma);

	//! Transposes the array.
	void Transpose();

    //! Computes \f$l_2\f$-norm.
    double Norm2() const;

	//! Computes \f$l_1\f$-norm.
    double Norm1() const;

	//! Computes \f$l_p\f$-norm.
    double Norm(double p) const;

    //! Computes the Hamming norm of a vector. Only implemented for Boolean type.
    double HammingNorm();

	//! Non-destructive element access.
	T Get(size_t i, size_t j) const;

	//! Overwrites data.
    void Set(std::shared_ptr<T> data);

    //! Access methods whose purpose is to provide a common interface with sparse matrix classes.
    void Set(size_t i, size_t j, T val) { this->operator()(i,j) = val; }

	//! Returns a column.
    CDenseVector<T> GetColumn(size_t j) const;

    //! Sets a column.
	void SetColumn(size_t j, const CDenseVector<T>& col);

	//! Returns a row.
    CDenseVector<T> GetRow(size_t i) const;

    //! Sets a row.
    void SetRow(size_t i, const CDenseVector<T>& row);

	//! Element access.
	T& operator()(size_t i, size_t j);

	//! Sums two arrays.
    CDenseArray<T> operator+(const CDenseArray<T>& array) const;

    //! Divides two arrays element-wise.
    CDenseArray<T> operator/(const CDenseArray<T>& array) const;

	//! Multiplies two arrays pointwise.
    CDenseArray<T> operator^(const CDenseArray<T>& array) const;

	//! Adds a scalar to all elements.
    CDenseArray<T> operator+(const T& scalar) const;

	//! Subtracts a scalar from all elements.
    CDenseArray<T> operator-(const T& scalar) const;

	//! Subtracts two arrays.
    CDenseArray<T> operator-(const CDenseArray<T>& array) const;

	//! Multiplies the array with a scalar.
    CDenseArray<T> operator*(const T& scalar) const;

    //! Divides the array by a scalar.
    CDenseArray<T> operator/(const T& scalar) const;

    //! Multiplies the object with an array from the right.
    //template<class Array> Array operator*(const Array& array) const;

	//! Multiplies the object with an array from the right.
    CDenseArray<T> operator*(const CDenseArray<T>& array) const;

	//! Multiplies the object with an array from the right.
    CDenseVector<T> operator*(const CDenseVector<T>& vector) const;

    //! Matrix-vector multiplication.
    template<u_int n> CVector<T,n> operator*(const CVector<T,n>& vector) const;

	//! Computes the standard inner product.
    static double InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y);

    //! Computes inner products between columns.
    static CDenseVector<T> ColumwiseInnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y);

    //! Computes the inner product w.r.t. to a kernel.
    static double InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y, CMercerKernel<T>& kernel);

	//! Computes tensor product of two matrices.
    static CDenseArray<T> KroneckerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y);

	//! Returns the transpose of a matrix.
    static CDenseArray<T> Transpose(const CDenseArray<T>& x);

	//! Access number of cols.
    size_t NRows() const { return m_nrows; }

	//! Access number of cols.
    size_t NCols() const { return m_ncols; }

	//! Returns number of elements.
	size_t NElems() const { return m_nrows*m_ncols; }

	//! Get pointer to the data.
    std::shared_ptr<T> Data() const { return m_data; }

	//! In-place scalar multiplication.
	void Scale(T scalar);

    //! Scale column.
    CDenseArray<T> ScaleColumns(const CDenseVector<T>& s);

    //! In-place addition of a scalar.
    void Add(const T& scalar);

	//! Sums up all elements.
	T Sum() const;

	//! Empirical mean of matrix entries.
	T Mean() const { return Sum()/NElems(); }

	//! Empirical variance of matrix entries.
    T Variance();

	//! Median of matrix entries;
    T Median();

    //! Median of absolute deviations.
    T MAD();

	//! Minimum of matrix entries;
	T Min() const;

	//! Maximum of matrix entries;
	T Max() const;

	//! Forms the absolute value of all matrix entries.
	void Abs();

    //! Writes matrix to a stream.
	template <class U> friend std::ostream& operator << (std::ostream& os, const CDenseArray<U>& x);

    //! Writes a boolean matrix to a file stream.
    friend std::ofstream& operator << (std::ofstream& os, const CDenseArray<bool>& x);

    //! Writes an integer matrix to a file stream.
    friend std::ofstream& operator << (std::ofstream& os, const CDenseArray<int>& x);

    //! Writes a float matrix to a file stream.
    friend std::ofstream& operator << (std::ofstream& os, const CDenseArray<float>& x);

    //! Writes a unsigned long integer matrix to a file stream.
    friend std::ofstream& operator << (std::ofstream& os, const CDenseArray<size_t>& x);

    //! Writes a double matrix to a file stream.
    friend std::ofstream& operator << (std::ofstream& os, const CDenseArray<double>& x);

	//! Reads matrix from a stream.
	template <class U> friend std::istream& operator >> (std::istream& is, CDenseArray<U>& x);

    //! Reads matrix from a file stream.
    template <class U> friend std::ifstream& operator >> (std::ifstream& is, CDenseArray<U>& x);

    //! Writes a matrix to a file.
    bool WriteToFile(const char* filename);

    //! Reads a matrix from file.
    bool ReadFromFile(const char* filename);

	//! Normalizes the matrix.
	bool Normalize();

	/*! \brief Computes the trace of a matrix.
	 *
	 * \details Works for non-square matrices. In this case, the sum is taken over the principal diagonal.
	 */
	T Trace() const;

	/*! \brief Computes the determinant of a matrix.
	 *
	 * \details The current implementation can only handle \f$2\times 2\f$- and \f$3\times 3\f$ matrices.
	 */
	T Determinant() const;

	//! Returns the number of bytes of the data type.
    size_t SizeOf() const { return sizeof(T); }

    //! Returns the numerical type.
    ETYPE GetType();

    //! Bilinear interpolation.
    double InterpolateBilinearly(double u, double v) const;

    //! Matrix inversion.
    bool Invert();

    //! Typecast operator
    template<typename Array> operator Array() {

        Array result(m_nrows,m_ncols);

        for(size_t i=0; i<result.NRows(); i++) {

            for(size_t j=0; j<result.NCols(); j++)
                result.Set(i,j,this->Get(i,j));

        }

        return result;

    }

protected:

	size_t m_nrows;				//!< number of rows
	size_t m_ncols;				//!< number of cols
	bool m_transpose;			//!< transpose flag
    std::shared_ptr<T> m_data;  //!< container that holds the array data

};

/*! \brief dense vector
 *
 *
 *
 */
template<class T>
class CDenseVector:public CDenseArray<T> {
public:

	//! Standard constructor.
	CDenseVector();

    /*! \brief Inherited constructor.
     *
     * \details Use this constructor to explicitly create row vectors.
     */
    CDenseVector(size_t nrows, size_t ncols);

    //! Constructs a column vector length \f$n\f$.
	CDenseVector(size_t n);

	//! Copy constructor.
	CDenseVector(const CDenseVector& vector);

    //! Constructs a column vector length \f$n\f$ from the given data.
    CDenseVector(size_t n, std::shared_ptr<T> data);

    //! Copy constructor.
    template<u_int n> CDenseVector(CVector<T,n>& x);

    //! Assignment operator.
    CDenseVector<T> operator=(const CDenseVector<T>& array);

    //! Deep copy.
    CDenseVector<T> Clone();

	//! Adds a scalar to a vector.
    CDenseVector<T> operator+(const T& scalar) const;

	//! Sums two vectors.
    CDenseVector<T> operator+(const CDenseVector<T>& vector) const;

    //! In-place addition of a vector.
    void Add(const CDenseVector<T>& vector);

	//! Subtracts two vectors.
    CDenseVector<T> operator-(const CDenseVector<T>& vector) const;

    //! Divides two vectors element-wise.
    CDenseVector<T> operator/(const CDenseVector<T>& vector) const;

	//! Multiplies the vector with a scalar.
    CDenseVector<T> operator*(const T& scalar) const;

	//! Element access.
	T& operator()(size_t i);

	//! Non-destructive element access.
	T Get(size_t i) const;

	//! Non-destructive element access.
    T Get(size_t i, size_t j) const { return CDenseArray<T>::Get(i,j); }

	//! Element access.
    T& operator()(size_t i, size_t j) { return CDenseArray<T>::operator ()(i,j); }

	//! Computes cross product of two 3-vectors.
    static CDenseVector<T> CrossProduct(const CDenseVector<T>& x, const CDenseVector<T>& y);

	//! Sorts the elements of the vector in ascending order.
	void Sort();

protected:

	using CDenseArray<T>::m_ncols;
	using CDenseArray<T>::m_nrows;
	using CDenseArray<T>::m_transpose;
	using CDenseArray<T>::m_data;

};

/*! \brief dense symmetric 2d matrix/array
 *
 *
 *
 */
template<class T>
class CDenseSymmetricArray {

public:

	//! Standard constructor.
	CDenseSymmetricArray();

	//! Constructor.
	CDenseSymmetricArray(size_t nrows);

	//! Copy constructor.
	CDenseSymmetricArray(const CDenseSymmetricArray& array);

	//! Destructor.
	~CDenseSymmetricArray();

	//! Print the array to standard ouput.
	void Print() const;

	//! \copydoc CDenseArray::Norm2()
	T Norm2() const;

	//! Element access.
	T& operator()(size_t i, size_t j);

	//! Non-destructive element access.
	T Get(size_t i, size_t j) const;

	//! Assignment operator.
	CDenseSymmetricArray<T> operator=(const CDenseSymmetricArray<T>& array);

	//! Sums two arrays.
	CDenseSymmetricArray<T> operator+(const CDenseSymmetricArray<T>& array) const;

	//! Subtracts two arrays.
	CDenseSymmetricArray<T> operator-(const CDenseSymmetricArray<T>& array) const;

	//! Multiplies the array with a scalar.
	CDenseSymmetricArray<T> operator*(const T& scalar) const;

	//! Multiplies the object with an array from the right.
	CDenseArray<T> operator*(const CDenseArray<T>& array) const;

	//! Access number of cols.
    size_t NRows() const { return m_nrows; }

	//! Access number of cols.
    size_t NCols() const { return m_nrows; }

	//! Get pointer to the data.
    T* Data() const { return m_data; }

	//! In-place scalar multiplication.
	void Scale(T scalar);

	//! Method stump (transpose does not do anything for symmetric matrices).
    void Transpose() {}

protected:

	size_t m_nrows;				//!< number of rows
	T* m_data;					//!< container that holds the array data

};


}

#endif /* DARRAY_H_ */
