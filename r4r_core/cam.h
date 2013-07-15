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

#ifndef R4RCAM_H_
#define R4RCAM_H_

#include <stdlib.h>
#include "trafo.h"
#include "types.h"

namespace R4R {

/*! \brief abstract camera class
 *
 *
 *
 */
class CAbstractCam {

public:

    /*! \brief Projects a point into the image plane.
     *
     * \param[in] x point in camera coordinates
     * \returns point in pixel coordinates
     *
     */
    virtual vec2 Project(const vec3& x) const = 0;
    virtual vec2f Project(const vec3f& x) const = 0;

    /*! \brief Projects a set of points into the image plane.
     *
     * \param[in] x points in camera coordinates
     * \returns points in pixel coordinates
     *
     */
    std::vector<vec2> Project(const std::vector<vec3>& x) const;
    std::vector<vec2f> Project(const std::vector<vec3f>& x) const;

    /*! \brief Computes the projection of a point into the image plane and its Jacobian.
     *
     * \param[in] x point in camera coordinates
     * \param[out] u point in pixel coordinates
     * \param[out] J Jacobian of the projection mapping in u
     *
     */
    virtual void Project(const vec3& x, vec2& u, mat& J) const = 0;
    virtual void Project(const vec3f& x, vec2f& u, matf& J) const = 0;

    /*! \brief Converts a pixel into a viewing direction.
     *
     * \param[in] u location w.r.t. the pixel coordinate system
     * \returns direction vector, normalized s.t. \f$z\f$-component equals \f$1\f$
     *
     */
    virtual vec3 Normalize(const vec2& u) const = 0;
    virtual vec3f Normalize(const vec2f& u) const = 0;

    /*! \brief Converts a set of pixels into viewing directions.
     *
     * \param[in] u locations w.r.t. the pixel coordinate system
     * \returns direction vectors
     *
     */
    std::vector<vec3> Normalize(const std::vector<vec2>& u) const;
    std::vector<vec3f> Normalize(const std::vector<vec2f>& u) const;

    /*! \brief Projects differential motion to optical flow vector.
     *
     * \param[in] x point in camera coordinates
     * \param[in] dx differential motion vector in 3d
     * \returns direction vectors
     *
     */
    virtual vec2 Flow(const vec3& x, const vec3& dx) const = 0;
    virtual vec2f Flow(const vec3f& x, const vec3f& dx) const = 0;

    //! Batch computation of optical flow.
    std::vector<vec2> Flow(const std::vector<vec3>& x, const std::vector<vec3>& dx) const;
    std::vector<vec2f> Flow(const std::vector<vec3f>& x, const std::vector<vec3f>& dx) const;

    //! Writes the camera parameters to a stream.
    virtual void Write(std::ostream& os) const = 0;

    //! Reads the camera parameters from a stream.
    virtual void Read(std::istream& is) = 0;

    //! Checks if two cameras are the same.
    virtual bool operator==(CAbstractCam& cam) = 0;

protected:

};


/*! \brief pinhole camera
 *
 *
 *
 */
class CPinholeCam:public CAbstractCam {

public:

    //! Constructor.
    CPinholeCam();

    //! Constructor.
    CPinholeCam(size_t w, size_t h);

    //! Constructor.
    CPinholeCam(double fu, double fv, double cu, double cv);

    //! \copydoc CAbstractCamera::Project(const vec3&) const
    vec2 Project(const vec3& x) const;
    vec2f Project(const vec3f& x) const;

    //! \copydoc CAbstractCamera::Project(const vec3&,vec2&,mat&) const
    void Project(const vec3& x, vec2& u, mat& J) const;
    void Project(const vec3f& x, vec2f& u, matf& J) const;

    //! \copydoc CAbstractCamera::Normalize()
    vec3 Normalize(const vec2& u) const;
    vec3f Normalize(const vec2f& u) const;

    //! \copydoc CAbstractCamera::Flow(const vec3&,const vec3&) const
    vec2 Flow(const vec3& x, const vec3& dx) const;
    vec2f Flow(const vec3f& x, const vec3f& dx) const;

    //! \copydoc CAbstractCamera::Write(std::ostream&) const
    void Write(std::ostream& os) const { os << *this; }

    //! Writes the camera parameters to a stream.
    friend std::ostream& operator << (std::ostream& os, const CPinholeCam& x);

    //! \copydoc CAbstractCamera::Read(std::istream&)
    void Read(std::istream& is) { is >> *this; }

    //! Reads the camera parameters from a stream.
    friend std::istream& operator >> (std::istream& is, CPinholeCam& x);

    //! Checks if two cameras are the same.
    bool operator==(CAbstractCam& cam);

    //! Access to image size.
    CVector<size_t,2> GetSize() {  return { m_size[0], m_size[1] }; }

    //! Projection matrix.
    mat GetProjectionMatrix();

private:

    size_t m_size[2];			//!< pixel size
    double m_f[2];				//!< focal length
    double m_c[2];				//!< principle point
    double m_alpha;				//!< skew coefficient
    double m_k[5];				//!< distortion coefficients

};

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
 *	of type CAbstractCam) and a pose (object of type CTransformation).
 *
 */
template<typename T=double>
class CView {

public:

    //! Constructor.
    CView(CAbstractCam& cam);

    //! Constructor.
    CView(CAbstractCam& cam, const CRigidMotion<T,3>& F);

    //! Copy constructor.
    CView(const CView<T>& view);

    //! Assignment operator.
    CView<T> operator=(const CView<T>& view);

    /*! \brief Projects a point into the image plane.
     *
     * \param[in] x point in world coordinates
     * \returns point in pixel coordinates
     *
     */
    CVector<T,2> Project(const CVector<T,3>& x);

    /*! \brief Projects a set of points into the image plane.
     *
     * \param[in] x vector of points in world coordinates
     * \returns points in pixel coordinates
     *
     */
    std::vector<CVector<T,2> > Project(const std::vector<CVector<T,3> >& x);

    /*! \brief Computes the projection of a point into the image plane and its Jacobian.
     *
     * \param[in] x point in world coordinates
     * \param[out] u point in pixel coordinates
     * \param[out] J Jacobian of the projection mapping in u
     *
     */
    void Project(const CVector<T,3>& x, CVector<T,2>& u, CDenseArray<T>& J);

    /*! \brief Converts a pixel into a viewing direction.
     *
     * \param[in] u location w.r.t. the pixel coordinate system
     * \returns direction vector, normalized s.t. \f$z\f$-component equals \f$1\f$ and transformed
     * to world coordinates
     *
     */
    CVector<T,3> Normalize(const CVector<T,2>& u);

    /*! \brief Converts a set of pixels into viewing directions.
     *
     * \param[in] u locations w.r.t. the pixel coordinate system
     * \returns direction vectors, normalized s.t. \f$z\f$-component equals \f$1\f$ and transformed
     * to world coordinates
     *
     */
    std::vector<CVector<T,3> > Normalize(const std::vector<CVector<T,2> >& u);

    //! Computes flow from differential motion in world coordinates.
    CVector<T,2> Flow(const CVector<T,3>& x, const CVector<T,3>& dx);

    //! Computes flow from differential motion in world coordinates.
    std::vector<CVector<T,2> > Flow(const std::vector<CVector<T,3> >& x,const std::vector<CVector<T,3> >& dx);

    /*! \brief Reads a intrisic/extrinsic parameters from a file.
     *
     * \param[in] filename file name
     *
     */
    bool OpenFromFile(const char* filename);

    /*! \brief Writes intrisic/extrinsic parameters to file.
     *
     * \param[in] filename file name
     *
     */
    bool SaveToFile(const char* filename);

    //! Writes the viewpoint parameters to a stream.
    template<typename U> friend std::ostream& operator << (std::ostream& os, const CView<U>& x);

    //! Access to the transformation.
    CRigidMotion<T,3>& GetTransformation() { return m_F; }

    //! Access to the inverse transformation.
    CRigidMotion<T,3>& GetInverseTransformation() { return m_Finv; }

    //! Access to the cam.
    CAbstractCam& GetCam() { return m_cam; }



protected:

    CAbstractCam& m_cam;               //!< intrinsic parameters
    CRigidMotion<T,3> m_F;             //!< extrinsic parameters
    CRigidMotion<T,3> m_Finv;          //!< inverse extrinsic parameters

};


}



#endif /* CAM_H_ */
