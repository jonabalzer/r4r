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

#ifndef R4RUTILS_H_
#define R4RUTILS_H_

#include "cam.h"
#include "types.h"

namespace R4R {

/*! \brief linear algebra utility functions
 *
 */
class CLinearAlgebra {

public:

	/*! \brief Computes the parameters of a Givens rotation in a numerical stable way.
	 *
	 * \param[in] a value above the element to eliminate
	 * \param[in] b value to eliminate
	 * \param[out] c cosine component of rotation matrix
	 * \param[out] s sine component of rotation matrix
	 * \param[out] r norm of (a,b)
 	 *
	 */
	template<class T> static void GivensRotation(T a, T b, T& c, T& s, T& r);

	//! Assembles a skew-symmetric matrix from input.
	static mat Hat(vec omega);

	/*! \brief Computes a rotation matrix from an axis-angle representation.
	 *
	 * \details Transcribed from VLFeat.
	 *
	 */
	static mat Rodrigues(vec omega);

	/*! \brief Computes a rotation matrix from an axis-angle representation.
	 *
	 * \details Transcribed from VLFeat.
	 *
	 */
	static mat Rodrigues(double o1, double o2, double o3);

	/*! Computes a rotation matrix from an axis-angle representation and its derivative.
	 *
	 * \details Transcribed from VLFeat.
	 *
	 */
	static void Rodrigues(vec omega, mat& R, mat& DRx, mat& DRy, mat& DRz);

	/*! \brief Implements a coordinate transformation of points.
	 *
	 * \param[in] F homogeneous transformation matrix
	 * \param[in] x non-homogeneous coordinate vector
	 * \returns point w.r.t. to new coordinate system
	 *
	 */
	static vec TransformPoint(mat& F, vec& x);

	/*! \brief Implements a coordinate transformation of directions.
	 *
	 * \param[in] F homogeneous transformation matrix
	 * \param[in] x non-homogeneous coordinate vector
	 * \returns direction w.r.t. to new coordinate system
	 *
	 */
	static vec TransformDirection(mat& F, vec& x);

	//! Computes a transformation matrix from translation and rotation parameters.
	static mat CreateTransformationMatrix(double tx, double ty, double tz, double ox, double oy, double oz);

	//! Computes a transformation matrix from translation and rotation parameters.
	static mat CreateInverseTransformationMatrix(double tx, double ty, double tz, double ox, double oy, double oz);

	//! Computes a transformation matrix from translation and rotation parameters.
    static mat CreateTransformationMatrix(vec params) { return CreateTransformationMatrix(params(0),params(1),params(2),params(3),params(4),params(5)); }

	//! Computes a transformation matrix from translation and rotation parameters.
	static mat CreateTransformationMatrix(mat& R, vec& t);

	/*! \brief Implements an inverse coordinate transformation of points.
	 *
	 * \param[in] F homogeneous transformation matrix
	 * \param[in] x non-homogeneous coordinate vector
	 * \returns point w.r.t. to new coordinate system
	 *
	 */
	static vec TransformPointBack(mat& F, vec& x);

	/*! \brief Implements an inverse coordinate transformation of directions.
	 *
	 * \param[in] F homogeneous transformation matrix
	 * \param[in] x non-homogeneous coordinate vector
	 * \returns direction w.r.t. to new coordinate system
	 *
	 */
	static vec TransformDirectionBack(mat& F, vec& x);

	//! Inverts a transformation.
	static mat InvertTransformation(mat& F);

	//! Extracts rotational part of transformation.
	static mat GetRotationMatrix(mat& F);

	//! Extracts translational part of transformation.
    static vec GetTranslationVector(mat& F);

    //! Projects a matrix onto SO(3) by Procrustes analysis.
	static mat ProjectToSO3(mat& R);

	/*! \brief Estimates the essential matrix by a DLTN technique.
	 *
	 * \param[in] corr pixel correspondences
	 * \param[in] cam intrinsic camera parameters
	 *
	 * \returns essential matrix
	 *
	 *
	 * \details As customary, the computed essential matrix describes the coordinate transformation from the first to the
	 * second frame (the former being assumed the world coordinate system).
	 *
	 */
	static mat CalibratedNPoint(const std::vector<std::pair<vec,vec> >& corr, const CCam& cam);

	/*! \brief Estimates a homography matrix \f$H:\mathrm{P}^2\to \mathrm{P}^2\f$ by a DLTN technique.
	 *
	 * \param[in] corr correspondences
	 * \returns homography matrix
	 *
	 * \details Make sure that input point sets are mean-free and normalized to improve quality of the result.
	 *
	 */
	static mat EstimateHomography(const std::vector<std::pair<vec,vec> >& corr);


	/*! \brief Decomposes an essential matrix into the only physically plausible frame.
	 *
	 * \param[in] E essential matrix
	 *
	 * \returns coordinate transformation from the first to the second frame
	 *
	 */
	static mat FactorEssentialMatrix(const mat& E);

	/*! \brief Extracts motion part from homography given intrinsic camera parameters.
	 *
	 * \param[in] H essential matrix
	 * \param[in] cam intrinsic camera parameters
	 *
	 * \returns coordinate transformation from the first to the second frame
	 *
	 */
	static mat ZhangFactorization(const mat& H, const CCam& cam);

};




}


#endif /* UTILSR4R_H_ */
