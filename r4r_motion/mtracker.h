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

#include "cam.h"
#include "stracker.h"
#include "lm.h"
#include "pcl.h"

namespace R4R {

class CMotionTrackerTracklet:public CSimpleTrackerTracklet {

    friend class CMotionTracker;

public:

    //! Deleted constructor.
    CMotionTrackerTracklet() = delete;

    //! Constructor.
    CMotionTrackerTracklet(size_t t0, const imfeature& x0, std::list<CInterestPoint<float,3> >::const_iterator it, size_t maxlength = 200):CSimpleTrackerTracklet(t0,x0,maxlength),m_pmap_point(it){}

private:

    std::list<CInterestPoint<float,3> >::const_iterator m_pmap_point;    //!< iterator directed to point in the map

};

/*! \brief motion tracker
 *
 *
*/
class CMotionTracker: public CSimpleTracker {

public:

    typedef CPointCloud<std::list,float,3> CSfMMap;

    //! Constructor.
    CMotionTracker(const CParameters* params, CPinholeCam<float>& cam);

    /*! \copybrief CTracker::Update(cv::Mat&,cv::Mat&)
     *
     *
     *
     */
    void Update(const std::vector<cv::Mat>& pyramid0, const std::vector<cv::Mat>& pyramid1);

    //! \copydoc CSimpleTracker::AddTracklets(cv::Mat&)
    void AddTracklets(const std::vector<cv::Mat>& pyramid);

    //! Returns a reference to the point cloud.
    const CSfMMap& GetMap() { return m_map; }

    //! Access to the motion.
    const C3dTrajectory<std::list,float>& GetMotion() { return m_motion; }

    //! Computes the last view.
    CView<float> GetLatestView();

private:

    CPinholeCam<float> m_cam;                    //!< camera
    C3dTrajectory<std::list,float> m_motion;     //!< motion
    CSfMMap m_map;                               //!< map

};

class CMagicSfM:public CLeastSquaresProblem<CCSRMatrix<float>,float> {

public:

	//! Constructor.
    CMagicSfM(CPinholeCam<float> cam, std::pair<std::vector<vec2f>,std::vector<vec2f> >& corri2i, std::pair<std::vector<vec3f>,std::vector<vec2f> >& corrs2i, CRigidMotion<float,3> F0inv);

	//! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
    void ComputeResidual(vecf& r) const;

	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&,const vec&)
    void ComputeResidualAndJacobian(vecf& r, CCSRMatrix<float>& J) const;

protected:

    CPinholeCam<float>& m_cam;											//!< intrinsic camera parameters
    std::pair<std::vector<vec2f>,std::vector<vec2f> >& m_corri2i;    	//!< image-to-image correspondences
    std::pair<std::vector<vec3f>,std::vector<vec2f> >& m_corrs2i;		//!< scene-to-image correspondences
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
