/*
 * factor.h
 *
 *  Created on: Jun 11, 2012
 *      Author: jbalzer
 */

#ifndef R4RFACTOR_H_
#define R4RFACTOR_H_

#include "darray.h"

namespace R4R {

/*! \brief LAPACK wrapper
 *
 */
template <class T>
class CMatrixFactorization {

public:

	/*! \brief Singular value decomposition.
	 *
	 * \param[in] A matrix to decompose
	 *
	 * \details
	 *
	 */
    static bool SVD(const CDenseArray<T>& A, CDenseArray<T>& U, CDenseArray<T>& S, CDenseArray<T>& Vt);

	//! In-place Cholesky decomposition of a symmetric matrix.
    static bool Cholesky(CDenseArray<T>& A);

	//! In-place matrix inversion by Cholesky decomposition.
    static bool InvertSymmetric(CDenseArray<T>& A);

	//! Rank of matrix.
    static size_t Rank(const CDenseArray<T>& A, T tol);

private:

};


}

#endif /* FACTOR_H_ */
