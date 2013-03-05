/*
 * mtracker.h
 *
 *  Created on: Jul 26, 2012
 *      Author: jbalzer
 */

#ifndef MTRACKER_H_
#define MTRACKER_H_

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
	CMotionTracker(CParameters params, CCam cam);

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

	CCam m_cam;						//!< camera
	std::list<vec> m_motion;		//!< motion


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
	vec Triangulate(vec o0, vec d0, vec o1, vec d1, double eps);


	/*! \brief Epipolar image correspondence search.
	 *
	 * \param[in] u0 pixel in the first image
	 * \param[in] ep epipole
	 * \param[in] d epipolar search direction in the second image
	 * \param[out] u1 pixel in the first image corresponding to \f$u_0\f$
	 * \param[out] lambda depth of the scene point corresponding to \f$u_0\f$ w.r.t. to projection center of first frame
	 * \param[in] img0 first image
	 * \param[in] img1 second image
	 * \param[in] hsize half-size of search window
	 *
	 */
	vec EpipolarSearch(vec& ub, vec& d, cv::Mat& img0, cv::Mat& img1, size_t hsize, size_t nsteps, double eps);

};



class CMagicSfM:public CLeastSquaresProblem<smat> {

public:

	//! Constructor.
	CMagicSfM(CCam cam, std::vector<std::pair<vec,vec> >& corri2i, std::vector<std::pair<vec,vec> >& corrs2i, mat F0inv);

	//! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
	virtual void ComputeResidual(vec& r);

	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&,const vec&)
	virtual void ComputeResidualAndJacobian(vec& r, smat& J);

	//! \copydoc CLeastSquaresProblem::ComputeDispersion(vec&)
	virtual vec ComputeDispersion(const vec& r);

protected:

	CCam m_cam;															//!< intrinsic camera parameters
	std::vector<std::pair<vec,vec> >& m_corri2i;						//!< image-to-image correspondences
	std::vector<std::pair<vec,vec> >& m_corrs2i;						//!< scene-to-image correspondences
	mat m_F0inv;														//!< transformation from first frame of image pair to world coordinates

};



/*! \brief energy for direct motion estimation
 *
 *
 * \details The template parameter is preselected to be of dense matrix type, as the Jacobian of this problem will be dense.
 *
 */
class CSfMTrackerUpdate:public CLeastSquaresProblem<mat> {

public:

	//! Constructor.
	CSfMTrackerUpdate(CCam cam, vector<CTracklet*> m_map, cv::Mat& img0, cv::Mat& img1, size_t hsize);

	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&)
	virtual void ComputeResidualAndJacobian(vec& r, mat& J);

	//! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
	virtual void ComputeResidual(vec& r);

private:

	CCam m_cam;												//!< intrinsic camera parameters
	std::vector<CTracklet*> m_map;							//!< pointers to tracks carrying a 3d point estimate
	cv::Mat& m_img0;										//!< frame 0
	cv::Mat& m_img1;										//!< frame 1
	size_t m_hsize;											//!< half of the LK search window

};

}

#endif /* MTRACKER_H_ */
