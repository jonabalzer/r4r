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

#ifndef R4RSTRACKER_H_
#define R4RSTRACKER_H_

#include "tracker.h"
#include "descriptor.h"


#include <opencv2/opencv.hpp>

namespace R4R {

// forward declaration
class CSimpleTracker;
class CMotionTracker;


// only for testing to switch between different implementations
typedef CCircularTracklet mytracklet;

/*! \brief tracklets for simple tracker
 *
 */
class CSimpleTrackerTracklet:public mytracklet {

    friend class CSimpleTracker;
    friend class CMotionTracker;

public:

    //! Standard constructor.
    CSimpleTrackerTracklet() = delete;

    //! Constructor.
    CSimpleTrackerTracklet(size_t t0, const imfeature& x0, size_t maxlength = 200):CTracklet(t0,x0,maxlength), m_reference_feature(x0) {}

private:

    imfeature m_reference_feature;        //! reference feature

};

/*! \brief simple LK tracker
 *
 *
* \details The following parameters have to be supplied by means of an object of CParameters:
 *
 * - SCALE number of scale levels to track
 * - FEATURE_THRESHOLD detection threshold
 * - TERMINATION_TIME number of frames to process in video stream
 * - REFRESH_RATE number of frames between redetections
 * - TRACKING_HSIZE half window size (for low-level motion estimation)
 * - DESCRIPTOR_HSIZE half-window size for descriptor computation
 * - MINIMAL_FEATURE_HDISTANCE half of the minimal distance that new features keep to existing ones
 * - MAX_HAMMING_DISTANCE threshold on hamming distance between initial and current BRIEF descriptor above which track
 * is dropped
 * - LK_PYRAMID_LEVEL number of scales to use in low-level LK motion estimation
 * - MAX_ITER maximum number of Gauss-Newton steps to execute in LK motion estimation
 * - ACCURACY time derivative of objective below which to break off Gauss-Newton iteration
 * - LAMBDA relative weight of the spatial image derivatives impact to the optical flow estimation
 *
 * Features stem from FAST detection with non-maximum suppression.  A weak minimal distance constraint is enforced
 * through integral images. A BRIEF descriptor is computed for every feature which enables a comparison with the
 * initial configuration of the corresponding tracklet. This provides a simple but efficient mechanism for occlusion
 * detection.
 *
 *
 */
class CSimpleTracker: public CSlidingWindowTracker {

public:

	//! Constructor.
    CSimpleTracker(const CParameters* params);

	/*! \copydoc CTracker::Init(cv::Mat& img)
	 *
	 * \details First, an image pyramid of height \f$SCALE+1\f$ is built. Features are extracted
	 * at the coarsest scale. They form the initial values of a set of tracklets that are attached
	 * to the tracker.
	 *
	 */
    void Init(const std::vector<cv::Mat>& pyramid);

	//! \copydoc CTracker::Update(std::vector<cv::Mat>&,std::vector<cv::Mat>&)
    void Update(const std::vector<cv::Mat>& pyramid0, const std::vector<cv::Mat>& pyramid1);

	//! \copydoc CTracker::Update(std::vector<cv::Mat>&)
    void UpdateDescriptors(const std::vector<cv::Mat>& pyramid);

	//! \copydoc CTracker::Clean(cv::Mat&,cv::Mat&)
    void Clean(const std::vector<cv::Mat>& pyramid0, const std::vector<cv::Mat>& pyramid1);

	//! \copydoc CTracker::AddTracklets(cv::Mat&)
    void AddTracklets(const std::vector<cv::Mat>& pyramid);

protected:

	cv::FastFeatureDetector m_detector;						//!< feature detector

};

} // end of namespace

#endif /* STRACKER_H_ */
