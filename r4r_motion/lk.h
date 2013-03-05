/*
 * lk.h
 *
 *  Created on: Aug 13, 2012
 *      Author: jbalzer
 */

#ifndef LK_H_
#define LK_H_



#include <opencv2/opencv.hpp>
#include "lm.h"

//#include <tbb/blocked_range.h>
//#define TBB_USE_DEBUG 1

namespace R4R {

class CLukasKanade:public CLeastSquaresProblem<mat>  {


/*! \brief Lukas-Kanade motion estimation functional
 *
 *
 *
 */
public:

	//! Constructor.
	CLukasKanade(const vec u0, const size_t hsize, const cv::Mat& img0, const cv::Mat& img1);

	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&)
	void ComputeResidualAndJacobian(vec& r, mat& J);

	//! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&, Matrix&, const vec&)
	void ComputeResidualAndJacobian(vec& r, mat& J, const vec& weights) {};

	//! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
	void ComputeResidual(vec& r) {};

private:

	const vec m_u0;												//! pixel in frame 0
	const size_t m_hsize;											//!< half of the LK search window
	const cv::Mat& m_img0;										//!< frame 0
	const cv::Mat& m_img1;										//!< frame 1

};
/*
class CLowLevelTracking {

public:

	//! Levenberg-Marquardt optimization of the LK functional.
	static vec LukasKanade(vec& u0, vec& t0, size_t hsize, cv::Mat& img0, cv::Mat& img1, size_t maxiter, double eps, double lambda);

	//! Multi-threaded Levenberg-Marquardt optimization of the LK functional.
	static void LukasKanadeMT(cv::Mat& img0, cv::Mat& img1, std::vector<vec>& points0, std::vector<vec>& points1, std::vector<uchar>& status, size_t hsize, size_t maxiter, double eps, double lambda);

};


class CLKTrackerInvoker {


public:

	//! Constructor.
	CLKTrackerInvoker(const cv::Mat& img0, const cv::Mat& img1, const vec* points0, vec* points1, uchar* status, size_t hsize, size_t maxiter, double eps, double lambda);

	//! Standard interface to function object.
	void operator()(const tbb::blocked_range<size_t>& range) const;

private:

	const cv::Mat* m_img0;
	const cv::Mat* m_img1;
	const vec* m_points0;
	vec* m_points1;
	uchar* m_status;
	size_t m_hsize;
	size_t m_maxiter;
	double m_eps;
	double m_lambda;

};*/


}



#endif /* LK_H_ */
