/*
 * ptracker.cpp
 *
 *  Created on: Jan 31, 2013
 *      Author: jbalzer
 */



#include "ptracker.h"

using namespace cv;

namespace R4R {



CPlaneTracker::CPlaneTracker(CParameters* params, CCam cam, Mat& img0):
		CSimpleTracker(params),
		m_cam(cam),
		m_img0(img0),
		m_motion(),
		m_texture_domain(),
		m_H(Mat::eye(3,3, CV_64FC1)) {

		m_motion.push_back(vec(6));

}


bool CPlaneTracker::Init(std::vector<cv::Mat>& pyramid) {

	CSimpleTracker::Init(pyramid);

	m_texture_domain = CTracker::GetManualBoundingBox(pyramid[0]);

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		vector<shared_ptr<CTracklet> > todel;

		// mark for deletion
		for(it=at(s).begin(); it!=at(s).end(); it++) {

			vec x = (*it)->GetLatestLocationAtNativeScale();
			Point2i xcv = Point2i(x.Get(0),x.Get(1));

			if(!m_texture_domain.contains(xcv)) {

				todel.push_back(move(*it));

			}

		}

		// delete and remove from list
		for(size_t i=0; i<todel.size(); i++) {

			todel[i].reset();
			at(s).remove(todel[i]);

		}

	}

	return 0;

}


bool CPlaneTracker::Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1) {

	// do LK tracking at every time instance
	CSimpleTracker::Update(pyramid0,pyramid1);

    size_t kfr = (size_t)m_params->GetIntParameter("KEYFRAME_RATE");

	if(m_global_t>0 && (m_global_t)%kfr==0) {

		// collect correspondence, maybe allocate space
		vector<Point2f> src, tgt;

		list<shared_ptr<CTracklet> >::iterator it;

		for(size_t s=0; s<size(); s++) {

			for(it=at(s).begin(); it!=at(s).end(); it++) {

				if((*it)->GetStatus() && (*it)->size()>(size_t)kfr) {

					vec u1 = (*it)->GetLatestLocationAtNativeScale();
					vec u0 = (*it)->GetPastLocationAtNativeScale((size_t)kfr);
					//vec u0 = (*it)->front()->GetLocationAtNativeScale();


					src.push_back(Point2f(u0.Get(0),u0.Get(1)));
					tgt.push_back(Point2f(u1.Get(0),u1.Get(1)));

				}

			}

		}

		// estimate homography
		//vector<bool> mask(src.size());
		Mat H = findHomography(src,tgt,CV_RANSAC,3);
		m_H = H*m_H;

		// find occlusion map between consecutive frames, set as global variable for clean routine
/*

		Mat img0w;

		warpPerspective(pyramid0[0],img0w,H,Size(pyramid0[0].cols,pyramid0[0].rows));

		Mat out = Mat::zeros(m_img0.rows,m_img0.cols,CV_8UC1);

		CIntegralImage<size_t> intimage = CIntegralImage<size_t>(out.cols,out.rows);

		for(int i=0; i<out.rows; i++) {

			for(int j=0; j<out.cols; j++) {

				size_t diff = (size_t)fabs((double)img0w.at<unsigned char>(i,j)-((double)pyramid1[0].at<unsigned char>(i,j)));
				intimage.AddDensityFast(j,i,diff);


			}

		}

		intimage.Compute();

		for(int i=0; i<out.rows; i++) {

			for(int j=0; j<out.cols; j++) {

				size_t vi = intimage.EvaluateFast(j,i,25,25);

				if(vi<20000 && (double)img0w.at<unsigned char>(i,j)!=0)
					out.at<unsigned char>(i,j) = 255;
				else
					out.at<unsigned char>(i,j) = 0;



			}

		}

		namedWindow("Test");
		imshow("Test",out);
		waitKey(0);
*/




	}

	return 0;

}

cv::Mat CPlaneTracker::PushforwardImage(cv::Mat& bg, cv::Mat& fg) {

	Mat fgw;

	assert(fg.channels()==3);

	warpPerspective(fg,fgw,m_H,Size(fg.cols,fg.rows));

	for(int i=0; i<fgw.rows; i++) {

		for(int j=0; j<fgw.cols; j++) {

			Vec3b fgwv = fgw.at<Vec3b>(i,j);

			if(fgwv.val[0]==0 && fgwv.val[1]==0 && fgwv.val[2]==0) {

				fgw.at<Vec3b>(i,j) = bg.at<Vec3b>(i,j);

			}


		}

	}



	return fgw;

}

cv::Mat CPlaneTracker::GenerateOcclusionMap(cv::Mat& img, double eps) {

	Mat img0w;

	warpPerspective(m_img0,img0w,m_H,Size(m_img0.cols,m_img0.rows));

	Mat out = Mat::zeros(m_img0.rows,m_img0.cols,CV_8UC1);

	// integral image
	CIntegralImage<size_t> intimage = CIntegralImage<size_t>(out.cols,out.rows);

	for(int i=0; i<out.rows; i++) {

		for(int j=0; j<out.cols; j++) {

			Vec3b imgv = img.at<Vec3b>(i,j);
			Vec3b img0v = img0w.at<Vec3b>(i,j);

			size_t diff = 0;

			for(size_t k=0; k<2; k++)
				diff += (size_t)fabs((double)imgv.val[k]-(double)img0v.val[k]);

			intimage.AddDensityFast(j,i,diff);
			/*if(diff<eps && (img0v.val[0]!=0 || img0v.val[1]!=0 || img0v.val[2]!=0))
				out.at<unsigned char>(i,j) = 255;
			else
				out.at<unsigned char>(i,j) = 0;*/

		}

	}

	intimage.Compute();

	CDenseArray<size_t> err = CDenseArray<size_t>(out.rows,out.cols);

	for(int i=0; i<out.rows; i++) {

		for(int j=0; j<out.cols; j++) {

			Vec3b img0v = img0w.at<Vec3b>(i,j);

			size_t vi = intimage.EvaluateFast(j,i,25,25);

			err(i,j) = vi;
			if(vi<eps && (img0v.val[0]!=0 || img0v.val[1]!=0 || img0v.val[2]!=0))
				out.at<unsigned char>(i,j) = 255;
			else
				out.at<unsigned char>(i,j) = 0;



		}

	}

    //if(m_global_t==141)
    //	err.SaveToFile("intimage.txt");

	return out;

}



}
