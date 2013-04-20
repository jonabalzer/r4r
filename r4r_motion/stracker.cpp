/*
 * stracker.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: jbalzer
 */

#include "stracker.h"
#include "basic.h"
#include "descspecial.h"
#include "feature.h"
#include "lk.h"
#include <fstream>

using namespace std;
using namespace cv;

namespace R4R {

CSimpleTracker::CSimpleTracker(CParameters* params):
	CTracker(params),
    m_detector(m_params->GetDoubleParameter("FEATURE_THRESHOLD"))
	{

	// generate sample points for tests performed in BRIEF descriptor
    CBRIEF::GenerateSamplePoints();

}

bool CSimpleTracker::Init(vector<Mat>& pyramid) {

	AddTracklets(pyramid);

	return 0;

}

bool CSimpleTracker::Update(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

	for(size_t s = 0; s<pyramid0.size(); s++) {

		// make a smooth copy of current pyramid level for descriptor computation
		Mat imsmooth;
        GaussianBlur(pyramid1[s],imsmooth,Size(0,0),m_params->GetDoubleParameter("GRAD_SMOOTH_SIGMA"));

		vector<Point2f> points0;
		vector<Point2f> points1;
		vector<uchar> status;
		vector<float> error;

		list<shared_ptr<CTracklet> >::iterator it;

		// get current feature locations
		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

				vec u0 = (*it)->GetLatestLocation();

				points0.push_back(Point2f(u0(0),u0(1)));
				points1.push_back(Point2f(u0(0),u0(1)));

			}

		}

		// check if there are features at that scale
		if(points0.size()>0) {

			// perform LK tracking
			calcOpticalFlowPyrLK(pyramid0[s],
								 pyramid1[s],
								 points0,
								 points1,
								 status,
								 error,
                                 Size(2*m_params->GetIntParameter("TRACKING_HSIZE")+1,2*m_params->GetIntParameter("TRACKING_HSIZE")+1),
                                 m_params->GetIntParameter("LK_PYRAMID_LEVEL"),
								 TermCriteria(TermCriteria::COUNT|TermCriteria::EPS,
                                              m_params->GetIntParameter("MAX_ITER"),
                                              m_params->GetDoubleParameter("ACCURACY")),
                                 m_params->GetDoubleParameter("LAMBDA"),
								 0);

			size_t counter = 0;

			// update tracklets
			for(it=at(s).begin(); it!=at(s).end(); it++) {

				if((*it)->GetStatus()) {

					// if status is ok and point in image, update feature
					if(status[counter] && points1[counter].x>0 && points1[counter].x<pyramid1[s].cols-1 && points1[counter].y>0 && points1[counter].y<pyramid1[s].rows-1) {

                        CFeature x(points1[counter].x,points1[counter].y,s);
						(*it)->Update(x);

					}
					else
						(*it)->SetStatus(false);

					counter++;

				}

			}

		}


	}


	m_global_t++;

	return 0;

}


void CSimpleTracker::Clean(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

	size_t no = 0;

	for(size_t s = 0; s<pyramid0.size(); s++) {

		// compute integral image for counting features
		CIntegralImage<size_t> cimg = ComputeFeatureDensity(pyramid0[0].cols/pow(2,s),pyramid0[0].rows/pow(2,s),s);
		cimg.Compute();

		list<shared_ptr<CTracklet> >::iterator it;

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			// if feature is still alive
			if((*it)->GetStatus()) {

				vec x = (*it)->GetLatestLocation();

				// check quality and distance criterion
                if((*it)->GetLatestState().GetQuality()>m_params->GetIntParameter("MAX_HAMMING_DISTANCE") ||
                   cimg.EvaluateFast(x(0),x(1),m_params->GetIntParameter("MINIMAL_FEATURE_HDISTANCE_CLEAN"),m_params->GetIntParameter("MINIMAL_FEATURE_HDISTANCE_CLEAN"))>1) {

					// delete last feature
					(*it)->pop_back();
					(*it)->SetStatus(false);

				} else
					no++;	// count number of valid tracks

			}

		}

	}

    m_n_active_tracks = no;

}

bool CSimpleTracker::AddTracklets(vector<Mat>& pyramid) {

	for(size_t s = 0; s<pyramid.size(); s++) {

		// make a smooth copy of current pyramid level for descriptor computation
		Mat imsmooth;
        GaussianBlur(pyramid[s],imsmooth,Size(0,0),m_params->GetDoubleParameter("GRAD_SMOOTH_SIGMA"));

		// compute integral image for counting features
		CIntegralImage<size_t> cimg = ComputeFeatureDensity(pyramid[s].cols,pyramid[s].rows,s);

		vector<KeyPoint> keypoints;
		m_detector.detect(pyramid[s],keypoints);

		if(keypoints.size()>0) {

			for(size_t i=0; i<keypoints.size(); i++)
				cimg.AddDensityFast(keypoints[i].pt.x,keypoints[i].pt.y,1);

			cimg.Compute();

            double hsize = m_params->GetIntParameter("MINIMAL_FEATURE_HDISTANCE_INIT");

			for(size_t i=0; i<keypoints.size(); i++) {

				if(cimg.EvaluateFast(keypoints[i].pt.x,keypoints[i].pt.y,hsize,hsize)==1) {

                    CFeature x(keypoints[i].pt.x,keypoints[i].pt.y,s);

					// create new tracklet with feature
					AddTracklet(x);

				}

			}

		}


	}

	return 0;

}

bool CSimpleTracker::UpdateDescriptors(std::vector<cv::Mat>& pyramid) {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		Mat imsmooth;
        GaussianBlur(pyramid[s],imsmooth,Size(0,0),m_params->GetDoubleParameter("GRAD_SMOOTH_SIGMA"));

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

                CFeature& x = (*it)->GetLatestState();
                vec u0 = x.GetLocation();

				// create new feature
				CRectangle<double> droi(u0.Get(0),
										u0.Get(1),
                                        m_params->GetIntParameter("DESCRIPTOR_HSIZE"),
                                        m_params->GetIntParameter("DESCRIPTOR_HSIZE"));

				// adjust region to scale
				droi.Scale(1.0/pow(2,s));

				// compute brief descriptor
                CBRIEF* briefdesc1 = new CBRIEF(droi);
                shared_ptr<CAbstractDescriptor> brief(briefdesc1);
				brief->Compute(imsmooth);
				//brief->Compute(pyramid[s]);
                x.AttachDescriptor("BRIEF",brief);

				// compute quality
                shared_ptr<CAbstractDescriptor> desc0 = (*it)->front().GetDescriptor("BRIEF");
				double quality = 1;

				if(desc0!=nullptr) {

                    shared_ptr<CBRIEF> briefdesc0 = static_pointer_cast<CBRIEF>(desc0);
					quality = briefdesc0->Distance(*briefdesc1);

				}

                x.SetQuality(quality);

                if(m_params->GetIntParameter("COMPUTE_ID")) {

					// compute identity descriptor
					CIdentityDescriptor* tempid = new CIdentityDescriptor(droi,
                                                                          (size_t)m_params->GetIntParameter("NORMALIZE_ID"),
                                                                          (size_t)m_params->GetIntParameter("DESCRIPTOR_HSIZE"));

					tempid->Compute(pyramid[s]);

                    shared_ptr<CAbstractDescriptor> id(tempid);
                    x.AttachDescriptor("ID",id);

				}

                if(m_params->GetIntParameter("COMPUTE_GRAD")) {

                    // compute normalized gradient field
                    CIdentityGradientDescriptor* temp = new CIdentityGradientDescriptor(droi,
                                                                                        m_params->GetDoubleParameter("ALPHA_GRAD_NORM"),
                                                                                        (size_t)m_params->GetIntParameter("NORMALIZE_GRAD"),
                                                                                        (size_t)m_params->GetIntParameter("DESCRIPTOR_HSIZE"));
                    temp->Compute(imsmooth);

                    shared_ptr<CAbstractDescriptor> idg(temp);
                    x.AttachDescriptor("GRAD",idg);

                }

                if(m_params->GetIntParameter("COMPUTE_HOG")) {

                    // compute HoG
                    CHistogramOfGradients* temp = new CHistogramOfGradients(droi);

                    temp->Compute(imsmooth);

                    shared_ptr<CAbstractDescriptor> pdesc(temp);
                    x.AttachDescriptor("HOG",pdesc);

                }

			}


		}

	}

	return 0;

}



}
