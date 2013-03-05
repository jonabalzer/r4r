/*
 * factor.h
 *
 *  Created on: Jun 11, 2012
 *      Author: jbalzer
 */

#ifndef FACTOR_H_
#define FACTOR_H_

#include "darray.h"

namespace R4R {

/*! \brief LAPACK wrapper
 *
 */
class CMatrixFactorization {

public:

	/*! \brief Singular value decomposition.
	 *
	 * \param[in] A matrix to decompose
	 *
	 * \details
	 *
	 */
	static bool SVD(const CDenseArray<double>& A, CDenseArray<double>& U, CDenseArray<double>& S, CDenseArray<double>& Vt);

	//! In-place Cholesky decomposition of a symmetric matrix.
	static bool Cholesky(CDenseArray<double>& A);

	//! In-place matrix inversion by Cholesky decomposition.
	static bool InvertSymmetric(CDenseArray<double>& A);

	//! Rank of matrix.
	static size_t Rank(const CDenseArray<double>& A, double tol);

private:

};


}

#endif /* FACTOR_H_ */
