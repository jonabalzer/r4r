/*
 * ptracker.h
 *
 *  Created on: Jan 31, 2013
 *      Author: jbalzer
 */

#ifndef PTRACKER_H_
#define PTRACKER_H_

#include "cam.h"
#include "kfilter.h"
#include "stracker.h"
#include "lm.h"

namespace R4R {

/*! \brief plane tracker
 *
 *
*/
class CPlaneTracker: public CSimpleTracker {

public:

	//! Constructor.
    CPlaneTracker(CParameters* params, CCam cam, cv::Mat& img0);

	/*! \copybrief CTracker::Init(std::vector<cv::Mat>&)
	 *
	 *
	 *
	 */
	bool Init(std::vector<cv::Mat>& pyramid);

	/*! \copybrief CTracker::Update(cv::Mat&,cv::Mat&)
	 *
	 *
	 *
	 */
	bool Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1);

	//! Cuts out the initial region of the image and warps according to the state of the motion tracker.
	cv::Mat PushforwardImage(cv::Mat& bg, cv::Mat& fg);

	//! Access to RoI.
	cv::Rect GetRoI() { return m_texture_domain; };

	//! Generates an occlusion map.
	cv::Mat GenerateOcclusionMap(cv::Mat& img, double eps);

private:

	CCam m_cam;						//!< camera
	cv::Mat m_img0;						//!< initial frame
	std::list<vec> m_motion;		//!< motion
	cv::Rect m_texture_domain;		//!< domain in the initial image replaced and pushed forward during tracking
	cv::Mat m_H;					//!< homography from the initial image to the current image

};

}


#endif /* PTRACKER_H_ */
