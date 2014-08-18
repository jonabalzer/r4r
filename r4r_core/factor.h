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
