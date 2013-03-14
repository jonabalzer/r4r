/*
 * darray.h
 *
 *  Created on: Apr 5, 2012
 *      Author: jbalzer
 */

#ifndef DARRAY_H_
#define DARRAY_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include "types.h"

namespace R4R {

static enum {  B1U = 1,
               C1U = 2, C1S = 3,
               S2U = 4, S2S = 5,
               I4S = 6, I4U = 7,
               F4S = 8,
               L8S = 9, L8U = 10,
               D8S = 11 } ETYPE;


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
	CDenseArray(size_t nrows, size_t ncols, T val = 0);

	//! Copy constructor.
	CDenseArray(const CDenseArray& array);

	//! Copy constructor.
	CDenseArray(size_t nrows, size_t ncols, T* data);

	//! Destructor.
	virtual ~CDenseArray();

	//! Creates an identity matrix.
	void Eye();

	//! Fills in ones.
	void Ones();

	//! Files the array with uniformly distributed random numbers between \f$0\f$ and \f$1\f$.
	void Rand();

	//! Transposes the array.
	void Transpose();

    //! Computes \f$l_2\f$-norm.
	T Norm2() const;

	//! Computes \f$l_1\f$-norm.
	T Norm1() const;

	//! Computes \f$l_p\f$-norm.
	T Norm(size_t p) const;

	//! Non-destructive element access.
	T Get(size_t i, size_t j) const;

	//! Overwrites data.
	void Set(T* data);

	//! Returns a column.
	CDenseVector<T> GetColumn(size_t j) const;

	//! Returns a column.
	void SetColumn(size_t j, const CDenseVector<T>& col);

	//! Returns a row.
	CDenseVector<T> GetRow(size_t i) const;

	//! Element access.
	T& operator()(size_t i, size_t j);

	//! Assignment operator.
	CDenseArray<T> operator=(const CDenseArray<T>& array);

	//! Sums two arrays.
    CDenseArray<T> operator+(const CDenseArray<T>& array) const;

	//! Multiplies two arrays pointwise.
    CDenseArray<T>& operator^(const CDenseArray<T>& array) const;

	//! Adds a scalar to all elements.
	CDenseArray<T> operator+(const T& scalar) const;

	//! Subtracts a scalar from all elements.
	CDenseArray<T> operator-(const T& scalar) const;

	//! Subtracts two arrays.
	CDenseArray<T> operator-(const CDenseArray<T>& array) const;

	//! Multiplies the array with a scalar.
	CDenseArray<T> operator*(const T& scalar) const;

	//! Multiplies the array with a scalar.
	CDenseArray<T> operator/(const T& scalar) const;

	//! Multiplies the object with an array from the right.
	CDenseArray<T> operator*(const CDenseArray<T>& array) const;

	//! Multiplies the object with an array from the right.
	CDenseVector<T> operator*(const CDenseVector<T>& vector) const;

	//! Computes the standard inner product.
	T static InnerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y);

	//! Computes tensor product of two matrices.
	CDenseArray<T> static KroneckerProduct(const CDenseArray<T>& x, const CDenseArray<T>& y);

	//! Returns the transpose of a matrix.
	CDenseArray<T> static Transpose(const CDenseArray<T>& x);

	//! Access number of cols.
	size_t NRows() const { return m_nrows; };

	//! Access number of cols.
	size_t NCols() const { return m_ncols; };

	//! Returns number of elements.
	size_t NElems() const { return m_nrows*m_ncols; }

	//! Get pointer to the data.
	T* Data() const { return m_data; };

	//! In-place scalar multiplication.
	void Scale(T scalar);

	//! Sums up all elements.
	T Sum() const;

	//! Empirical mean of matrix entries.
	T Mean() const { return Sum()/NElems(); }

	//! Empirical variance of matrix entries.
	T Variance() const;

	//! Median of matrix entries;
	T Median() const;

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
	size_t SizeOf() { return sizeof(T); };

    //! Computes the Hamming norm of a vector. Only implemented for Boolean type.
    double HammingNorm();

protected:


	size_t m_nrows;				//!< number of rows
	size_t m_ncols;				//!< number of cols
	bool m_transpose;			//!< transpose flag
	T* m_data;					//!< container that holds the array data


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

	//! Inherited constructor.
	CDenseVector(size_t nrows, size_t ncols);

	//! Parametized constructor.
	CDenseVector(size_t n);

	/*! \brief Constructor.
	 *
	 * \details The constructed vector is a column vector by default. Use CDenseVector::CDenseVector(size_t n, bool row)
	 * to get a row vector.
	 *
	 * \param[in] n number of elements
	 * \param[in] val initial values
	 */
	//CDenseVector(size_t n, T val = 0);

	//! Constructor.
	//CDenseVector(size_t n, bool row);

	//! Copy constructor.
	CDenseVector(const CDenseVector& vector);

	//! Copy constructor.
	CDenseVector(size_t n, T* data, bool row = false);

	//! Adds a scalar to a vector.
	CDenseVector<T> operator+(const T& scalar) const;

	//! Sums two vectors.
	CDenseVector<T> operator+(const CDenseVector<T>& vector) const;

	//! Subtracts two vectors.
	CDenseVector<T> operator-(const CDenseVector<T>& vector) const;

	//! Multiplies the vector with a scalar.
	CDenseVector<T> operator*(const T& scalar) const;

	//! Element access.
	T& operator()(size_t i);

	//! Non-destructive element access.
	T Get(size_t i) const;

	//! Non-destructive element access.
	T Get(size_t i, size_t j) const { return CDenseArray<T>::Get(i,j); };

	//! Element access.
	T& operator()(size_t i, size_t j) { return CDenseArray<T>::operator ()(i,j); };

	//! Computes cross product of two 3-vectors.
	CDenseVector<T> static CrossProduct(const CDenseVector<T>& x, const CDenseVector<T>& y);

	//! Sorts the elements of the vector in ascending order.
	void Sort();

protected:

	using CDenseArray<T>::m_ncols;
	using CDenseArray<T>::m_nrows;
	using CDenseArray<T>::m_transpose;
	using CDenseArray<T>::m_data;

};


/*! \brief convenience class for creating vectors of length \f$n\f$
 *
 *
 *
 */
template<class T, int n>
class CNVector: public CDenseVector<T> {

public:

    //! Constructor.
    CNVector():CDenseVector<T>(n) {};


};

/*! \brief convenience class for creating vectors of length \f$3\f$
 *
 * \todo What about operators in base class>
 *
 */
template<class T>
class CNVector<T,3>:public CDenseVector<T> {

public:

    // constructor
    CNVector():CDenseVector<T>(3) {};

    // constructor
    CNVector(T x, T y, T z):CDenseVector<T>(3) {};

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
	size_t NRows() const { return m_nrows; };

	//! Access number of cols.
	size_t NCols() const { return m_nrows; };

	//! Get pointer to the data.
	T* Data() const { return m_data; };

	//! In-place scalar multiplication.
	void Scale(T scalar);

	//! Method stump (transpose does not do anything for symmetric matrices).
	void Transpose() {};

protected:

	size_t m_nrows;				//!< number of rows
	T* m_data;					//!< container that holds the array data

};


}

#endif /* DARRAY_H_ */
