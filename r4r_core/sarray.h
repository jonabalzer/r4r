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

#ifndef R4RSARRAY_H_
#define R4RSARRAY_H_

#include <map>
#include <vector>
#include <stdlib.h>
#include "darray.h"

namespace R4R {

/*!
 * \brief matrix-market CSR triple
 */
template<typename T,typename U = size_t>
class CCSRTriple {

public:

    //! Constructor.
    CCSRTriple(U i, U j, T v):m_i(i),m_j(j),m_v(v){}

    //! Lexicographic ordering.
    bool operator<(const CCSRTriple& x) const;

    //! Tests for equal location.
    bool operator==(const CCSRTriple& x) const { return m_i==x.m_i && m_j==x.m_j; }

    //! Tests for equal location.
    bool operator!=(const CCSRTriple& x) const { return m_i!=x.m_i || m_j!=x.m_j; }

    //! Access to row index.
    const U& i() const { return m_i; }

    //! Access to column index.
    const U& j() const { return m_j; }

    //! Access to value.
    const T& v() const { return m_v; }

private:

    U m_i;
    U m_j;
    T m_v;

};

/*!
 * \brief matrix-market CSR triple
 */
template<typename T,typename U = size_t>
class CCSCTriple {

public:

    //! Constructor.
    CCSCTriple(U i, U j, T v):m_i(i),m_j(j),m_v(v){}

    //! Initializer list constructor.
    //CCSCTriple(std::initializer_list<U> list);

    //! Lexicographic ordering.
    bool operator<(const CCSCTriple& x) const;

    //! Tests for equal location.
    bool operator==(const CCSCTriple& x) const { return m_i==x.m_i && m_j==x.m_j; }

    //! Tests for equal location.
    bool operator!=(const CCSCTriple& x) const { return m_i!=x.m_i || m_j!=x.m_j; }

    //! Access to row index.
    const U& i() const { return m_i; }

    //! Access to column index.
    const U& j() const { return m_j; }

    //! Access to value.
    const T& v() const { return m_v; }

private:

    U m_i;
    U m_j;
    T m_v;

};


// forward declaration
template<class T,typename U> class CCSCMatrix;

/*! \brief sparse matrix in compressed-row format
 *
 */
template<typename T,typename U = size_t>
class CCSRMatrix {

public:

    //! Deleted standard constructor.
    CCSRMatrix() = delete;

    //! Constructor.
    CCSRMatrix(size_t m, size_t n);

    /*! \brief Constructor which takes MatrixMarket triples as input.
     *
     * The triples can contain duplicate entries. A heap sort and implicit summation
     * is performed internally in the constructor.
     *
     */
    CCSRMatrix(size_t m, size_t n, std::vector<CCSRTriple<T,U> >& data);

    /*! \brief Constructor for external assembly.
     *
     * This facilitates fast assembly, e.g., in the case of least-squares
     * problems, where one iterates linearly over the abscissae, which
     * corresponds to a row, and for each abscissa, all the variables which
     * it affects (the columns). No random access is needed. Make sure to
     * call Verify() to see if internally everything is in order.
     *
     */
    CCSRMatrix(std::shared_ptr<std::vector<U> >& rowptr, std::shared_ptr<std::vector<U> >& cols, std::shared_ptr<std::vector<T> >& vals);

    //! Verifies the structure.
    bool Verify() const;

    //! Access number of cols.
    size_t NRows() const { return m_nrows; }

    //! Access number of cols.
    size_t NCols() const { return m_ncols; }

    //! In-place scalar multiplication.
    void Scale(T scalar);

    //! Erases the matrix and replaces it with the identity.
    void Eye();

    //! Counts the number of non-zero entries.
    size_t NNz() const { return m_vals->size(); }

    //! Multiplies the object with a dense array from the right.
    template<class Matrix> Matrix operator*(const Matrix& array) const;

    //! Writes matrix to a stream.
    template<typename V,typename W> friend std::ostream& operator << (std::ostream& os, const CCSRMatrix<V,W>& x);

    //! Transposition transform a CSR into a CSC matrix.
    static CCSRMatrix<T,U> Transpose(const CCSRMatrix<T,U>& x);

protected:

    size_t m_nrows;                                         //!< number of rows
    size_t m_ncols;                                         //!< number of cols
    bool m_transpose;                                       //!< transposition flag
    std::shared_ptr<std::vector<U> > m_rowptr;              //!< indicates the beginning of rows in #m_val and #m_cols
    std::shared_ptr<std::vector<U> > m_cols;                //!< col index
    std::shared_ptr<std::vector<T> > m_vals;                //!< value

};


template<typename T,typename U = size_t>
class CSymmetricCSRMatrix {

public:

    //! Constructor.
    CSymmetricCSRMatrix(size_t s);

    //! Squares a CSC matrix.
    CSymmetricCSRMatrix(const CCSCMatrix<T,U>& x);

    //! Access number of cols.
    size_t NRows() const { return m_size; }

    //! Access number of cols.
    size_t NCols() const { return m_size; }

    //! In-place scalar multiplication.
    void Scale(T scalar);

    //! Counts the number of non-zero entries.
    size_t NNz() const;

    //! Multiplies the object with a dense array from the right.
    template<class Matrix> Matrix operator*(const Matrix& array) const;

    //! Writes matrix to a stream.
    template<typename V,typename W> friend std::ostream& operator << (std::ostream& os, const CSymmetricCSRMatrix<V,W>& x);

private:

    size_t m_size;                                          //!< number of rows and cols
    std::shared_ptr<std::vector<U> > m_rowptr;              //!< indicates the beginning of rows in #m_val and #m_cols
    std::shared_ptr<std::vector<U> > m_cols;                //!< col index
    std::shared_ptr<std::vector<T> > m_vals;                //!< value

};

template<typename T,typename U = size_t>
class CCSCMatrix {

    friend class CSymmetricCSRMatrix<T,U>;

public:

    //! Deleted standard constructor.
    CCSCMatrix() = delete;

    //! Constructor.
    CCSCMatrix(size_t m, size_t n);

    /*! \brief Constructor which takes MatrixMarket triples as input.
     *
     * The triples can contain duplicate entries. A heap sort and implicit summation
     * is performed internally in the constructor.
     *
     */
    CCSCMatrix(size_t m, size_t n, std::vector<CCSCTriple<T,U> >& data);

    /*! \brief Constructor for external assembly.
     */
    CCSCMatrix(std::shared_ptr<std::vector<U> >& colptr, std::shared_ptr<std::vector<U> >& rows, std::shared_ptr<std::vector<T> >& vals);

    //! Erases the matrix and replaces it with the identity.
    void Eye();

    //! Access number of cols.
    size_t NRows() const { return m_nrows; }

    //! Access number of cols.
    size_t NCols() const { return m_ncols; }

    //! In-place scalar multiplication.
    void Scale(T scalar);

    //! Counts the number of non-zero entries.
    size_t NNz() const { return m_vals->size(); }

    /*! Multiplies the transpose of the current object with an array from the right.
     *
     * CAVEAT: This implicitly form the transpose of the current object!
     *
     */
    template<class Matrix> Matrix operator*(const Matrix& array) const;

    //! Writes matrix to a stream.
    template<typename V,typename W> friend std::ostream& operator << (std::ostream& os, const CCSCMatrix<V,W>& x);

    /*! \brief Squares a matrix in an efficient way.
     *
     * This comes in handy for forming normal equations. The entries of
     * the squares are obtained as all mutual inner product between the
     * columns, so the representation has no influence on speed. The CSC
     * matrix can be kept for fast multiplication of the right-hands side
     * of a linear system, while the squared matrix in CSR  format is particularly
     * fit for any Krylov subspace method which only needs matrix-vector-products
     * (without transposition), e.g. the standard CG method. Problem is however
     * that assembly is usually done row-wise! So assembly is followed by
     * transposition (inefficient), than squaring. Maybe it is faster to do
     * construct CSC matrix from triplets and a heap.
     *
     *
     */
    CSymmetricCSRMatrix<T,U> Square() const;

    //! Transposition transform a CSR into a CSC matrix.
    static CCSCMatrix<T,U> Transpose(const CCSCMatrix<T,U>& x);

private:

    size_t m_nrows;                                         //!< number of rows
    size_t m_ncols;                                         //!< number of cols
    bool m_transpose;                                       //!< transposition flag
    std::shared_ptr<std::vector<U> > m_colptr;              //!< indicates the beginning of columns in #m_val and #m_rows
    std::shared_ptr<std::vector<U> > m_rows;                //!< row index
    std::shared_ptr<std::vector<T> > m_vals;                //!< value

};

// forward declarations
template <class T> class CSparseDiagonalArray;
template <class T> class CSparseLowerTriangularArray;
template <class T> class CSparseUpperTriangularArray;

/*! \brief sparse 2d matrix/array
 *
 *
 */
template<class T>
class CSparseArray {

	friend class CSparseDiagonalArray<T>;
	friend class CSparseLowerTriangularArray<T>;
	friend class CSparseUpperTriangularArray<T>;

    typedef std::map<size_t,std::map<size_t,T> > spdata;

public:

	//! Constructor
	CSparseArray();

	//! Constructor
	CSparseArray(size_t nrows, size_t ncols);

    //! Copy constructor
    CSparseArray(const CSparseArray<T>& x);

    //! Copy constructor.
    CSparseArray(size_t nrows, size_t ncols, std::shared_ptr<spdata> data);

    //! Assignment operator for shallow copies.
    CSparseArray<T>& operator =(const CSparseArray<T>& x);

    //! Makes a deep copy.
    CSparseArray<T> Clone() const;

    //! Access to the data container.
    std::shared_ptr<spdata> Data() { return m_data; }

    //! Creates an identiy matrix.
    void Eye();

    //! Resizes the array.
    void Resize(size_t nrows, size_t ncols) { m_nrows = nrows; m_ncols = ncols; }

    //! Concatenates the array with another.
    void Concatenate(const CSparseArray& array, bool direction);

	//! Writes matrix to a stream.
	template<class U> friend std::ostream& operator << (std::ostream& os, CSparseArray<U>& x);

	//! Computes the standard L2 norm.
	T Norm2();

	//! Transposes the array in place.
    void Transpose();

	//! Returns a transposed copy of the array.
	CSparseArray<T> static Transpose(const CSparseArray<T>& array);

	//! Efficiently computes the square \f$A^{\top}A\f$ of a sparse matrix \f$A\f$.
    CSparseArray<T> static Square(const CSparseArray<T>& array);

	//! Non-destructive element access.
    T Get(size_t i, size_t j) const;

	//! Returns a row of the matrix.
    std::map<size_t,T> GetRow(size_t i) const;

	/*! \brief Element access.
	 *
	 *	\details Opposed to CSparseArray::Get(size_t i, size_t j) and CSparseArray::Set(size_t i, size_t j, T v),
	 *	this routine causes fill-in, as accessing a non-existing element of the map container creates
	 *	this a new element and returns a reference to it.
	 *
	 * \param[in] i row index
	 * \param[in] j column index
	 * \return reference to the element
	 *
	 */
	T& operator()(size_t i, size_t j);

	//! Sets a matrix element to a given value.
	void Set(size_t i, size_t j, T v);

	//! Multiplies the array with a scalar.
	CSparseArray<T> operator*(const T& scalar);

	//! Sums two arrays.
	CSparseArray<T> operator+(const CSparseArray<T>& array);

	//! Subtracts two arrays.
	CSparseArray<T> operator-(const CSparseArray<T>& array);

	//! Multiplies the object with a dense array from the right.
    template<class Matrix> Matrix operator*(const Matrix& array) const;

    //! Computes the standard inner product.
    T static InnerProduct(const CSparseArray<T>& x, const CSparseArray<T>& y);

	//! Access number of cols.
    size_t NRows() const { return m_nrows; }

	//! Access number of cols.
    size_t NCols() const { return m_ncols; }

	//! In-place scalar multiplication.
	void Scale(T scalar);

    //! In-place scaling of row.
    void ScaleRow(size_t i, T scalar);

	//! Counts the number of non-zero entries.
	size_t Nonzeros();

	//! Checks whether the matrix is symmetric.
	bool Symmetric();

	//! Deletes a row.
	void DeleteRow(size_t i);

	//! Deletes a column.
	void DeleteColumn(size_t j);

    /*! Converts the matrix into compressed sparse row format suitable for many standard sparse solvers.
     *
     * TODO: Return R4R type.
     *
     */
	void GetCSR(std::vector<size_t>& nz, std::vector<size_t>& j, std::vector<T>& v, bool ibase);

	//! Converts the matrix into (row,column,value) format.
	void GetCOO(std::vector<size_t>& i, std::vector<size_t>& j, std::vector<T>& v, bool ibase);

	//! Saves the matrix in matrix market format.
    bool WriteToFile(const char* filename);

	//! Sets random entries in the matrix to random values.
	void Rand(size_t nnz);

    //! Typecast operator
    template<typename Array> operator Array() {

        Array result(m_nrows,m_ncols);

        for(size_t i=0; i<result.NRows(); i++) {

            for(size_t j=0; j<result.NCols(); j++)
                result(i,j) = this->Get(i,j);

        }

        return result;

    }

protected:

    size_t m_nrows;                      //!< number of rows
    size_t m_ncols;                      //!< number of cols
    bool m_transpose;                    //!< transpose flag
    std::shared_ptr<spdata> m_data;      //!< container that holds the array data

};

typedef CSparseArray<double> smat;
typedef CSparseArray<float> smatf;

/*! \brief sparse 2d matrix/array
 *
 *
 *
 */
template<class T>
class CSparseBandedArray {

public:

	//! Constructor
	CSparseBandedArray();

	//! Constructor
	CSparseBandedArray(size_t nrows, size_t ncols);

    //! Writes matrix to a stream.
    template<class U> friend std::ostream& operator << (std::ostream& os, const CSparseBandedArray<U>& x);

	//! Non-destructive element access.
	virtual T Get(size_t i, size_t j);

	//! Returns a range of bands.
	std::map<int,std::map<size_t,T> > GetBands(int lower, int upper);

	//! Sets a matrix element to a given value.
	virtual void Set(size_t i, size_t j, T v);

	/*! \brief Element access.
	 *
	 *
	 * \details Opposed to CSparseBandedArray::Get(size_t i, size_t j) and CSparseBandedArray::Set(size_t i, size_t j, T v),
	 * this routine causes fill-in, as accessing a non-existing element of the map container creates this a new element
	 * this and returns a reference to it.
	 *
	 * \param[in] i row index
	 * \param[in] j column index
	 *
	 * \return reference to the element
	 *
	 *
	 */
	virtual T& operator()(size_t i, size_t j);

	//! Access number of cols.
    size_t NRows() const { return m_nrows; }

	//! Access number of cols.
    size_t NCols() const { return m_ncols; }

	//! Computes the standard L2 norm.
	T Norm2();

	//! Transposes the array in place.
	virtual void Transpose();

	//! Returns a transposed copy of the array.
	CSparseBandedArray<T> static Transpose(const CSparseBandedArray<T>& array);

	/*! \brief Returns true if the array is transposed.
	 *
	 *
	 * \details No, transpose is an internal state flag only. Remove this routine for security reasons.
	 */
    bool IsTransposed() { return m_transpose; }

	//! In-place scalar multiplication.
	void Scale(T scalar);

	/*! \brief Add a scalar in-place.
	 *
	 * \details Do this only for the diagonal, otherwise the band structure will be lost.
	 *
	 */
	void AddDiagonal(const T scalar);

	//! Scale the diagonal in-place.
	void ScaleDiagonal(T scalar);

	//! Sums two arrays.
	CSparseBandedArray<T> operator+(const CSparseBandedArray<T>& array);

	//! Subtracts two arrays.
	CSparseBandedArray<T> operator-(const CSparseBandedArray<T>& array);

	//! Multiplies the object with a dense array from the right.
	CDenseArray<T> operator*(const CDenseArray<T>& array);

	//! Multiplies the object with a sparse array from the right.
	CSparseBandedArray<T> operator*(CSparseBandedArray<T>& array);

    //! Solves linear system associated with this matrix. FIXME: Implement this.
    virtual void Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const {}

	//! Deletes an element.
	void Delete(size_t i, size_t j);

	//! Deletes a row.
	void DeleteRow(size_t i);

	//! Deletes a column.
	void DeleteCol(size_t i);

	//! Counts the number of non-zeros.
	size_t Nonzeros();

	//! Returns the number of bands.
    size_t NoBands() { return m_data.size(); }

	//! Saves the matrix in matrix market format.
	bool SaveToFile(const char* filename);

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

	size_t m_nrows;														//!< number of rows
	size_t m_ncols;														//!< number of cols
	bool m_transpose;													//!< transpose flag
	std::map<int,std::map<size_t,T> > m_data;							//!< container that holds the array data

	//! Computes the band index.
    int Band(size_t i, size_t j) const { return (int)j - (int)i; }

	//! Computes the location in the band.
    size_t BandIndex(size_t i, size_t j) const { 	return j<=i ? j : i; }

	//! Computes the row from band/diagonal index representation.
    size_t Row(int b, size_t d) const { return b<=0 ? size_t((int)d-b) : d; }

	//! Computes the row from band/diagonal index representation.
    size_t Col(int b, size_t d) const { return b<=0 ? d : size_t((int)d+b); }

};

typedef CSparseBandedArray<double> sbmat;

/*! \brief sparse square diagonal matrix
 *
 *
 *
 */
template<class T>
class CSparseDiagonalArray: public CSparseBandedArray<T> {

public:

	//! Constructor.
	CSparseDiagonalArray();

	//! Constructor.
	CSparseDiagonalArray(size_t nrows, size_t ncols);

    //! Non-destructive element access.
    T Get(size_t i, size_t j);

	//! Non-destructive element access.
    T Get(size_t i) { return CSparseBandedArray<T>::Get(i,i); }

	//! \copydoc CSparseBandedArray::Set(size_t,size_t,T)
	void Set(size_t i, size_t j, T v);

	//! Sets a matrix element to a given value.
    void Set(size_t i, T v) { Set(i,i,v); }

	//! \copydoc CSparseBandedArray::operator()(size_t,size_t)
	T& operator()(size_t i, size_t j);

	//! \copydoc CSparseBandedArray::Solve(CDenseArray<T>&,const CDenseArray<T>&)
    void Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const;

	//! Method stump (transposition does nothing to square diagonal matrices).
    virtual void Transpose() {}

	//! Computes the inverse in-place.
	void Invert();

protected:

	using CSparseBandedArray<T>::m_data;
	using CSparseBandedArray<T>::m_nrows;
	using CSparseBandedArray<T>::m_ncols;

};

/*! \brief sparse upper-triangular matrix
 *
 * \todo Maybe one triangular banded class is enough + possibility to transpose?
 * \todo Check whether banded storage scheme is really faster for back-substitution.
 *
 */
template<class T>
class CSparseUpperTriangularArray: public CSparseBandedArray<T> {

public:

	//! Constructor.
	CSparseUpperTriangularArray();

	//! Constructor.
	CSparseUpperTriangularArray(size_t nrows, size_t ncols);

	//! \copydoc CSparseBandedArray::Get(size_t,size_t)
	T Get(size_t i, size_t j);

	//! \copydoc CSparseBandedArray::Set(size_t,size_t,T)
	void Set(size_t i, size_t j, T v);

	//! \copydoc CSparseBandedArray::operator()(size_t,size_t)
	T& operator()(size_t i, size_t j);

	//! \copydoc CSparseBandedArray::Solve(CDenseArray<T>&,const CDenseArray<T>&)
    void Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const;

	//! Method stump (no transposition allowed for triangular matrices).
    virtual void Transpose() {}

protected:

	using CSparseBandedArray<T>::m_data;
	using CSparseBandedArray<T>::m_nrows;
	using CSparseBandedArray<T>::m_ncols;
	using CSparseBandedArray<T>::m_transpose;

};


/*! \brief sparse lower-triangular matrix
 *
 *
 *
 */
template<class T>
class CSparseLowerTriangularArray: public CSparseBandedArray<T> {

public:

	//! Constructor.
	CSparseLowerTriangularArray();

	//! Constructor.
	CSparseLowerTriangularArray(size_t nrows, size_t ncols);

	//! \copydoc CSparseBandedArray::Get(size_t,size_t)
	T Get(size_t i, size_t j);

	//! \copydoc CSparseBandedArray::Set(size_t,size_t,T)
	void Set(size_t i, size_t j, T v);

	//! \copydoc CSparseBandedArray::operator()(size_t,size_t)
	T& operator()(size_t i, size_t j);

	//! \copydoc CSparseBandedArray::Solve(CDenseArray<T>&,const CDenseArray<T>&)
    void Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const;

    //! Method stump (no transposition allowed for triangular matrices).
    virtual void Transpose() {}

protected:

	using CSparseBandedArray<T>::m_data;
	using CSparseBandedArray<T>::m_nrows;
	using CSparseBandedArray<T>::m_ncols;
	using CSparseBandedArray<T>::m_transpose;


};

}

#endif /* SARRAY_H_ */
