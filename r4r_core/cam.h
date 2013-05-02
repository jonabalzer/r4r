//////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////////

#ifndef R4RCAM_H_
#define R4RCAM_H_

#include <stdlib.h>
#include "trafo.h"


namespace R4R {



/*! \brief camera model
 *
 *
 *
 */
class CCam {

public:

	size_t m_size[2];			//!< pixel size
	double m_f[2];				//!< focal length
	double m_c[2];				//!< principle point
	double m_alpha;				//!< skew coefficient
	double m_k[5];				//!< distortion coefficients
	mat m_F;					//!< frame world-to-cam
	mat m_Finv;					//!< inverse cam-to-world frame

	//! Constructor.
	CCam();

	//! Constructor.
	CCam(size_t w, size_t h);

	/*! \brief Projects a point into the image plane.
	 *
	 * \param[in] x point in world coordinates
	 * \returns point in pixel coordinates
	 *
	 */
	vec Project(vec x);

	/*! \brief Computes the projection of a point into the image plane and its Jacobian.
	 *
	 * \param[in] x point in world coordinates
	 * \param[out] u point in pixel coordinates
	 * \param[out] J Jacobian of the projection mapping in u
	 *
	 */
	void Project(vec x, vec& u, mat& J);

	/*! \brief Computes the projection of a point into the image plane and its Jacobian.
	 *
	 * \param[in] xc point in camera coordinates
	 * \param[out] u point in pixel coordinates
	 * \param[out] J Jacobian of the projection mapping in u
	 *
	 */
	void ProjectLocal(vec xc, vec& u, mat& J);

	/*! \brief Projects a point into the image plane without using the extrinsics.
	 *
	 * \param[in] xc point in camera coordinates
	 * \returns point in pixel coordinates
	 *
	 */
	vec ProjectLocal(vec xc);

	/*! \brief Projects a direction vector into the image plane without using the extrinsics.
	 *
	 * \param[in] d direction in camera coordinates
	 * \returns direction in camera coordinates
	 *
	 */
	vec ProjectDirectionLocal(vec dc);

	/*! \brief Converts a pixel into a viewing direction (in world coordinates).
	 *
	 * \param[in] u location w.r.t. the pixel coordinate system
	 * \returns direction vector, normalized s.t. z-component equals 1
	 *
	 */
	vec Normalize(vec u);

	/*! \brief Converts a pixel into a viewing direction (in camera coordinates).
	 *
	 * \param[in] u location w.r.t. the pixel coordinate system
	 * \returns direction vector, normalized s.t. z-component equals 1
	 *
	 */
	vec NormalizeLocal(const vec& u) const;

	/*! \brief Reads camera parameters from a file.
	 *
	 * \param[in] filename file name
	 *
	 */
	bool OpenFromFile(const char* filename);

	//! Writes camera parameters to file.
	bool SaveToFile(const char* filename);

	//! Writes the camera parameters to a stream.
	friend std::ostream& operator << (std::ostream& os, const CCam& x);

	//! Returns projection center in world coordinates.
	vec GetOrigin();

	//! Assembles projection matrix.
	mat GetProjectionMatrix() const;

	//! Assembles inverse projection matrix.
	mat GetInverseProjectionMatrix() const;

protected:


};


/*! \brief viewpoint
 *
 *	\details A viewpoint encompasses both, an (intrinsic) projection model (i.e., a reference to an object
 *	of type CCam) and a pose (object of type CTransformation).
 *
 */
class CView {

public:



protected:



};


}



#endif /* CAM_H_ */
