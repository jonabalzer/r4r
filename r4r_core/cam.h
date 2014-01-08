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

    //! Constructor.
    CPinholeCam(size_t w, size_t h, double fu, double fv, double cu, double cv);

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
    CVector<size_t,2> GetSize() const {  return { m_size[0], m_size[1] }; }

    //! Projection matrix.
    mat GetProjectionMatrix() const;

private:

    size_t m_size[2];			//!< pixel size
    double m_f[2];				//!< focal length
    double m_c[2];				//!< principle point
    double m_alpha;				//!< skew coefficient
    double m_k[5];				//!< distortion coefficients

};

/*! \brief viewpoint
 *
 *	\details A viewpoint encompasses both, an (intrinsic) projection model (i.e., a reference to an object
 *	of type CAbstractCam) and a pose (object of type CRigidMotion).
 *
 */
template<typename T=double>
class CView {

public:

    //! Deleted standard constructor. One always needs a reference to intrinsic camera model.
    CView() = delete;

    //! Constructor.
    CView(CAbstractCam& cam, u_int index = -1);

    //! Constructor.
    CView(CAbstractCam& cam, const CRigidMotion<T,3>& F, u_int index = -1);

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
    CVector<T,2> Project(const CVector<T,3>& x) const;

    /*! \brief Projects a set of points into the image plane.
     *
     * \param[in] x vector of points in world coordinates
     * \returns points in pixel coordinates
     *
     */
    std::vector<CVector<T,2> > Project(const std::vector<CVector<T,3> >& x) const;

    /*! \brief Computes the projection of a point into the image plane and its Jacobian.
     *
     * \param[in] x point in world coordinates
     * \param[out] u point in pixel coordinates
     * \param[out] J Jacobian of the projection mapping in u
     *
     */
    void Project(const CVector<T,3>& x, CVector<T,2>& u, CDenseArray<T>& J) const;

    /*! \brief Converts a pixel into a viewing direction.
     *
     * \param[in] u location w.r.t. the pixel coordinate system
     * \returns direction vector, normalized s.t. \f$z\f$-component equals \f$1\f$ and transformed
     * to world coordinates
     *
     */
    CVector<T,3> Normalize(const CVector<T,2>& u) const;

    /*! \brief Converts a set of pixels into viewing directions.
     *
     * \param[in] u locations w.r.t. the pixel coordinate system
     * \returns direction vectors, normalized s.t. \f$z\f$-component equals \f$1\f$ and transformed
     * to world coordinates
     *
     */
    std::vector<CVector<T,3> > Normalize(const std::vector<CVector<T,2> >& u) const;

    //! Computes flow from differential motion in world coordinates.
    CVector<T,2> Flow(const CVector<T,3>& x, const CVector<T,3>& dx) const;

    //! Computes flow from differential motion in world coordinates.
    std::vector<CVector<T,2> > Flow(const std::vector<CVector<T,3> >& x,const std::vector<CVector<T,3> >& dx) const;

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
     * The world-to-cam transformation is stored only. Upon reading, the inverse
     * can be re-computed.
     *
     */
    bool SaveToFile(const char* filename);

    //! Writes the viewpoint parameters to a stream.
    template<typename U> friend std::ostream& operator << (std::ostream& os, const CView<U>& x);

    //! Access to the transformation.
    const CRigidMotion<T,3>& GetTransformation() const { return m_F; }

    //! Access to the inverse transformation.
    const CRigidMotion<T,3>& GetInverseTransformation() const { return m_Finv; }

    //! Set transformation.
    void SetTransformation(const CRigidMotion<T,3>& F) { m_F = F; m_Finv = F; m_Finv.Invert(); }

    //! Set transformation.
    void SetInverseTransformation(const CRigidMotion<T,3>& Finv) { m_Finv = Finv; m_F = Finv; m_F.Invert(); }

    //! Access to the cam.
    const CAbstractCam& GetCam() const { return m_cam; }

    //! Gets the location of the projection center in world coordinates.
    CVector<T,3> GetLocation() const;

    //! Gets the principal axis of the projection device in world coordinates.
    CVector<T,3> GetPrincipalAxis() const;

    //! Return index of the view.
    int GetIndex() const { return m_index; }

    /*! Translates the origin of the view.
     *
     * \param[in] t translation vector in world coordinates
     *
     */
    void Translate(const CVector<T,3>& t);

    /*! Translates the origin of the view.
     *
     * \param[in] t translation vector in camera coordinates
     *
     */
    void DifferentialTranslate(const CVector<T,3>& t);

    /*! Orbits the camera around a point.
     *
     * \param[in] center center of orbit in world coordinates
     * \param[in] axis rotation axis in world coordinates
     */
    void Orbit(const CVector<T,3>& center, const CVector<T,3>& axis);

    //! Rotates the view so that it looks at a given point.
    void LookAt(CVector<T,3> x);

protected:

    CAbstractCam& m_cam;               //!< intrinsic parameters
    CRigidMotion<T,3> m_F;             //!< extrinsic parameters
    CRigidMotion<T,3> m_Finv;          //!< inverse extrinsic parameters
    int m_index;                       //!< view counter

};

/*! functor for comparing two views */
template <typename T=double>
class CIndexViewComparator {

    //! Compares two views based on their index.
    bool operator()(const CView<T>& x, const CView<T>& y) { return x.GetIndex() <= y.GetIndex(); }

};

/*! functor for comparing two views */
template <typename T=double>
class CDistanceViewComparator {

public:

    //! Constructor.
    CDistanceViewComparator(const CVector<T,3>& x): m_x(x) {}

    //! Rates two views based on the distance to a given point.
    bool operator()(const CView<T>& x, const CView<T>& y);

protected:

    CVector<T,3> m_x;  //! point in space

};

/*! functor for comparing two views */
template <typename T=double>
class CAngleViewComparator: public CDistanceViewComparator<T> {

public:

    //! Constructor.
    CAngleViewComparator(const CVector<T,3>& x, const CVector<T,3>& n): CDistanceViewComparator<T>::CDistanceViewComparator(x), m_n(n) { m_n.Normalize(); }

    //! Rates two views based on the angle of their principal axes with the normal at that point.
    bool operator()(const CView<T>& x, const CView<T>& y);

    /*! Computes the key for comparison.
     */
    T GetKey(const CView<T>& x);

private:

    CVector<T,3> m_n;  //! normal

    using CDistanceViewComparator<T>::m_x;

};

/*! functor for comparing two pairs of views */
template <typename T=double>
class CDistanceViewPairComparator {

public:

    //! Rates a pair based on the distance between its two views.
    bool operator()(const std::pair<CView<T>,CView<T> > & x, const std::pair<CView<T>,CView<T> >& y);

};

template <typename T=double>
class CAngleViewPairComparator {

public:

    //! Rates a pair based on the angle between their principal axes.
    bool operator()(const std::pair<CView<T>,CView<T> > & x, const std::pair<CView<T>,CView<T> >& y);

};

}



#endif /* CAM_H_ */
