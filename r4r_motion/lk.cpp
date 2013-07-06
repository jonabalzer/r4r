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

#include "lk.h"
//#include <tbb/tbb.h>
#include "interp.h"

using namespace std;
using namespace cv;
//using namespace tbb;

namespace R4R {

CLukasKanade::CLukasKanade(const vec u0, const size_t hsize, const cv::Mat& img0, const cv::Mat& img1):
	CLeastSquaresProblem((2*hsize+1)*(2*hsize+1),2),
	m_u0(u0),
	m_hsize(hsize),
	m_img0(img0),
	m_img1(img1) {

}

void CLukasKanade::ComputeResidualAndJacobian(vec& r, mat& J) {

	size_t w = 2*m_hsize+1;
	vec u1 = m_u0 + m_model;

	for(int i=-(int)m_hsize; i<=(int)m_hsize; i++) {

		for(int j=-(int)m_hsize; j<=(int)m_hsize; j++) {

			// compute row index
			int row = w*(i+(int)m_hsize)+(int)m_hsize+j;

			// interpolated intensities and gradients
			double I0, I1, I1u, I1v;
			I0 = CImageInterpolation::Bilinear(m_img0,m_u0.Get(0)+i,m_u0.Get(1)+j);
			I1 = CImageInterpolation::Bilinear(m_img1,u1(0)+i,u1(1)+j);
			I1u = CImageInterpolation::Gradient(m_img1,u1(0)+i,u1(1)+j,true);
			I1v = CImageInterpolation::Gradient(m_img1,u1(0)+i,u1(1)+j,false);

			// set residual
			r(row) = I1 - I0;

			// Jacobian
			J(row,0) = I1u;
			J(row,1) = I1v;

		}

	}

}

/*vec CLowLevelTracking::LukasKanade(vec& u0, vec& t0, size_t hsize, cv::Mat& img0, cv::Mat& img1, size_t maxiter, double eps, double lambda) {

	CLukasKanade problem(u0,hsize,img0,img1);

	vec& model = problem.Get();
	model = t0;

	mat M(0,0);
	CPreconditioner<mat,vec,double> precond = CPreconditioner<mat,vec,double>(M);
	CIterativeSolver<mat,vec,double> solver = CIterativeSolver<mat,vec,double>(precond,2,0,true);
	CLevenbergMarquardt<mat> lms(problem,solver,lambda);

	lms.Iterate(maxiter,eps,true);

	return u0 + model;

}

void CLowLevelTracking::LukasKanadeMT(cv::Mat& img0, cv::Mat& img1, vector<vec>& points0, vector<vec>& points1, vector<uchar>& status, size_t hsize, size_t maxiter, double eps, double lambda) {

	size_t n = points0.size();

	if(n==0)
		return;

	assert(points1.size()==n && status.size()==n);

	//parallel_for(blocked_range<size_t>(0,n), CLKTrackerInvoker(img0,img1,&points0[0],&points1[0],&status[0],hsize,maxiter,eps,lambda));


	vec* pts1 = &points1[0];
	uchar* st = &status[0];

	parallel_for( blocked_range<size_t>(0,n),[=](const blocked_range<size_t>& range) {

		for(size_t i=range.begin(); i!=range.end(); ++i) {

			// create problem instance
			CLukasKanade problem(points0[i],hsize,img0,img1);

			// init model
			vec& model = problem.Get();
			model = points1[i] - points0[i];

			// set up solver
			mat M(0,0);
			CPreconditioner<mat,vec,double> precond = CPreconditioner<mat,vec,double>(M);
			CIterativeSolver<mat,vec,double> solver = CIterativeSolver<mat,vec,double>(precond,2,0,true);
			CLevenbergMarquardt<mat> lms(problem,solver,lambda);

			// iterate
			lms.Iterate(maxiter,eps,true);

			vec u1 = points0[i] + model;

			// check whether points is out of bounds
			if(u1(0)>=5 && u1(1)<=img1.cols-6 && u1(1)>=5 && u1(1)<=img1.rows-6) {

				pts1[i] = u1;
				st[i] = 1;

			}
			else
				st[i] = 0;


		}


	});

}

CLKTrackerInvoker::CLKTrackerInvoker(const cv::Mat& img0, const cv::Mat& img1, const vec* points0, vec* points1, uchar* status, size_t hsize, size_t maxiter, double eps, double lambda):
		m_img0(&img0),
		m_img1(&img1),
		m_points0(points0),
		m_points1(points1),
		m_status(status),
		m_hsize(hsize),
		m_maxiter(maxiter),
		m_eps(eps),
		m_lambda(lambda) {

}

void CLKTrackerInvoker::operator()(const blocked_range<size_t>& range) const {

	// create local variables to speed up multithreading
	const Mat& img0 = *m_img0;
	const Mat& img1 = *m_img0;
	const vec* points0 = m_points0;
	vec* points1 = m_points1;
	uchar* status = m_status;

	for(size_t i=range.begin(); i!=range.end(); ++i) {

		// create problem instance
		CLukasKanade problem(points0[i],m_hsize,img0,img1);

		// init model
		vec& model = problem.Get();
		model = points1[i] - points0[i];

		// set up solver
		mat M(0,0);
		CPreconditioner<mat,vec,double> precond = CPreconditioner<mat,vec,double>(M);
		CIterativeSolver<mat,vec,double> solver = CIterativeSolver<mat,vec,double>(precond,2,0,true);
		CLevenbergMarquardt<mat> lms(problem,solver,m_lambda);

		// iterate
		vec r = lms.Iterate(m_maxiter,m_eps,true);

		vec u1 = points0[i] + model;

		// check whether points is out of bounds
		if(u1(0)>=5 && u1(1)<=img1.cols-6 && u1(1)>=5 && u1(1)<=img1.rows-6) {

			points1[i] = points1[i]-points0[i];

			cout << points0[i] << endl;
			status[i] = 1;
		}
		else
			status[i] = 0;

	}


}*/


}
