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

#ifndef R4RMTRACKER_H_
#define R4RMTRACKER_H_

#define COMPUTE_ID 0

#include "cam.h"
#include "kfilter.h"
#include "stracker.h"
#include "lm.h"

namespace R4R {

/*! \brief motion tracker
 *
 *
*/
class CMotionTracker: public CSimpleTracker {

public:

    //! Constructor.
    CMotionTracker(CParameters* params, CPinholeCam& cam);

    /*! \copybrief CTracker::Update(cv::Mat&,cv::Mat&)
     *
     *
     *
     */
    bool Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1);

    /*! \copybrief CTracker::Update(std::vector<cv::Mat>&)
     *
     * \param[in] pyramid color image pyramid
     *
     * Color is attached to the scene points as a descriptor and can be exported to the .PLY file.
     *
     */
    bool UpdateDescriptors(std::vector<cv::Mat>& pyramid);

    //! Computes color of map points.
    bool ColorMap(cv::Mat& img);

    //! Returns point cloud.
    std::map<CTracklet*,vec> GetMap();

    //! Exports the trajectory to file.
    bool SaveMotion(const char* filename);

    //! Exports the point cloud to file.
    bool SaveMap(const char* filename);

private:

    CPinholeCam m_cam;						//!< camera
    std::list<vec> m_motion;                //!< motion


    /*! \brief Triangulates the depth of a scene point from correspondence between images from fully calibrated cameras.
     *
     * \param[in] o0 location of first projection center
     * \param[in] d0 direction of first viewing ray
     * \param[in] o1 location of second projection center
     * \param[in] d1 direction of first viewing ray
     *
     * \details All vectors are represented in world coordinates.
     *
     */
    vec3f Triangulate(vec3f o0, vec3f d0, vec3f o1, vec3f d1, float eps);

    //! \brief Epipolar image correspondence search.
    vec2f EpipolarSearch(vec2f& ub, vec2f& d, cv::Mat& img0, cv::Mat& img1, size_t hsize, size_t nsteps, float eps);

};

class CMagicSfM:public CLeastSquaresProblem<smat,double> {

public:

	//! Constructor.
    CMagicSfM(CPinholeCam cam, std::vector<std::pair<vec2f,vec2f> >& corri2i, std::vector<std::pair<vec3f,vec2f> >& corrs2i, CRigidMotion<float,3> F0inv);

	//! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
    void ComputeResidual(vec& r);

	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&,const vec&)
    void ComputeResidualAndJacobian(vec& r, smat& J);

	//! \copydoc CLeastSquaresProblem::ComputeDispersion(vec&)
    vec ComputeDispersion(const vec& r);

protected:

    CPinholeCam& m_cam;													//!< intrinsic camera parameters
    std::vector<std::pair<vec2f,vec2f> >& m_corri2i;    				//!< image-to-image correspondences
    std::vector<std::pair<vec3f,vec2f> >& m_corrs2i;					//!< scene-to-image correspondences
    CRigidMotion<float,3> m_F0inv;										//!< transformation from first frame of image pair to world coordinates

};

///*! \brief energy for direct motion estimation
// *
// *
// * \details The template parameter is preselected to be of dense matrix type, as the Jacobian of this problem will be dense.
// *
// */
//class CSfMTrackerUpdate:public CLeastSquaresProblem<mat> {

//public:

//	//! Constructor.
//	CSfMTrackerUpdate(CCam cam, vector<CTracklet*> m_map, cv::Mat& img0, cv::Mat& img1, size_t hsize);

//	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&)
//	virtual void ComputeResidualAndJacobian(vec& r, mat& J);

//	//! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
//	virtual void ComputeResidual(vec& r);

//private:

//	CCam m_cam;												//!< intrinsic camera parameters
//	std::vector<CTracklet*> m_map;							//!< pointers to tracks carrying a 3d point estimate
//	cv::Mat& m_img0;										//!< frame 0
//	cv::Mat& m_img1;										//!< frame 1
//	size_t m_hsize;											//!< half of the LK search window

//};

} // end of namespace

#endif /* MTRACKER_H_ */
