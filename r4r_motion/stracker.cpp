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

#include "stracker.h"

#include <fstream>
#include <algorithm>

#include "descspecial.h"

using namespace std;
using namespace cv;

namespace R4R {

CSimpleTracker::CSimpleTracker(const CParameters *params):
	CTracker(params),
    m_detector(m_params->GetDoubleParameter("FEATURE_THRESHOLD"))
	{

	// generate sample points for tests performed in BRIEF descriptor
    CBRIEF::GenerateSamplePoints();

}

void CSimpleTracker::Init(const std::vector<Mat>& pyramid) {

    // detect and add tracklets
    AddTracklets(pyramid);

}

void CSimpleTracker::Update(const vector<Mat>& pyramid0, const vector<Mat>& pyramid1) {

    // sort features according to scale, this setup allows for scale changes
    vector<vector<Point2f> > points0(pyramid0.size());
    vector<vector<Point2f> > points1(pyramid0.size());
    vector<vector<shared_ptr<mytracklet> > > tracklets(pyramid0.size());

    list<shared_ptr<mytracklet > >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        // only consider live tracks
        if((*it)->GetStatus()) {

            const imfeature& f0 = (*it)->GetLatestState();
            const vec2f& u0 = f0.GetLocation();
            u_int scale = u_int(f0.GetScale());

            // make sure we can track at that level
            if(scale<pyramid0.size() && scale<pyramid1.size()) {

                points0[scale].push_back(Point2f(u0.Get(0),u0.Get(1)));
                points1[scale].push_back(Point2f(u0.Get(0),u0.Get(1)));
                tracklets[scale].push_back(*it);                            // keep order in which tracklets are added

            }

        }

    }

    for(size_t s=0; s<pyramid0.size(); s++) {

        // error/status flags
        vector<uchar> status;
		vector<float> error;

		// check if there are features at that scale
        if(points0[s].size()>0) {

			// perform LK tracking
			calcOpticalFlowPyrLK(pyramid0[s],
								 pyramid1[s],
                                 points0[s],
                                 points1[s],
								 status,
								 error,
                                 Size(2*m_params->GetIntParameter("TRACKING_HSIZE")+1,2*m_params->GetIntParameter("TRACKING_HSIZE")+1),
                                 m_params->GetIntParameter("LK_PYRAMID_LEVEL"),
								 TermCriteria(TermCriteria::COUNT|TermCriteria::EPS,
                                              m_params->GetIntParameter("MAX_ITER"),
                                              m_params->GetDoubleParameter("ACCURACY")),
                                 m_params->GetDoubleParameter("LAMBDA"),
                                 0);

            // update tracklets
            for(size_t i=0; i<tracklets[s].size(); i++) {

                if(status[i] && points1[s][i].x>0 && points1[s][i].x<pyramid1[s].cols-1 && points1[s][i].y>0 && points1[s][i].y<pyramid1[s].rows-1) {

                    vec2f loc = { points1[s][i].x, points1[s][i].y };
                    imfeature x(loc,s,0);
                    tracklets[s][i]->Update(x);

                }
                else
                    tracklets[s][i]->SetStatus(false);

            }

        }

    }

	m_global_t++;

}

void CSimpleTracker::Clean(const vector<Mat>& pyramid0, const vector<Mat>& pyramid1) {

    list<shared_ptr<mytracklet> >::iterator it;

    // now go through all tracklets
    for(it=m_data.begin(); it!=m_data.end(); it++) {

        const imfeature& f = (*it)->GetLatestState();

        // if feature is still alive  but its quality is too low in terms of the
        //   distance between its descriptor and the reference, then kill it
        if((*it)->GetStatus() && f.GetQuality()>m_params->GetIntParameter("MAX_HAMMING_DISTANCE"))
            (*it)->SetStatus(false);

	}

}

void CSimpleTracker::AddTracklets(const vector<Mat>& pyramid) {

    // initialize pyramid of integral images
    vector<CIntImage<size_t> > imgs;
    for(u_int s=0; s<pyramid.size(); s++)
        imgs.push_back(CIntImage<size_t>(pyramid[s].cols,pyramid[s].rows));

    // collect and count feature we already have per scale
    vector<size_t> n = ComputeFeatureDensity(imgs);
    size_t active = std::accumulate(n.begin(),n.end(),0);

    // only do something if we have too little tracks
    if(active<=size_t(m_params->GetIntParameter("MIN_NO_FEATURES"))) {

        int hsize = m_params->GetIntParameter("MINIMAL_FEATURE_HDISTANCE_INIT");
        vec2f hsize2 = { float(hsize), float(hsize) };

        for(u_int s = 0; s<pyramid.size(); s++) {

            // how many to add per scale
            int ntoadd = m_params->GetIntParameter("MAX_NO_FEATURES") - (int)n[s];

            // if we have enough, don't do anything
            if(ntoadd>0) {

                // detect
                vector<KeyPoint> keypoints;
                m_detector.detect(pyramid[s],keypoints);

                // only do something if there were detections at scale s
                if(keypoints.size()>0) {

                    // keep the best but no more than the number to add
                    KeyPointsFilter::retainBest(keypoints,ntoadd);

                    // add potential
                    for(size_t i=0; i<keypoints.size(); i++)
                        imgs[s].AddMass(keypoints[i].pt.y,keypoints[i].pt.x,1);

                    // compute the integral image
                    imgs[s].Compute();

                    for(size_t i=0; i<keypoints.size(); i++) {

                        vec2f loc = { float(keypoints[i].pt.x), float(keypoints[i].pt.y) };

                        // only features that are separated from all other candidates and existing features are accepted
                        if(imgs[s].EvaluateApproximately(loc,hsize2)==1) {

                            // create feature
                            imfeature x(loc,s,0);

                            // create new tracklet with feature, FIXME: get size restriction from parameters
                            CSimpleTrackerTracklet* tracklet = new CSimpleTrackerTracklet(m_global_t,x,m_params->GetIntParameter("BUFFER_LENGTH"));
                            this->AddTracklet(tracklet);

                            // also set the reference feature
                            tracklet->m_reference_feature = x;

                        }

                    }

                }

            }

        }

    }
    else {  // only compute integral images

        for(u_int s = 0; s<pyramid.size(); s++)
            imgs[s].Compute();

    }

    // since we already have the integral images computed, we might as well
    //  check whether two tracks got too close to each other
    list<shared_ptr<mytracklet> >::iterator it;

    int hsize = m_params->GetIntParameter("MINIMAL_FEATURE_HDISTANCE_CLEAN");
    vec2f hsize2 = { float(hsize), float(hsize) };

    // now go through all tracklets
    for(it=m_data.begin(); it!=m_data.end(); it++) {

        const imfeature& f = (*it)->GetLatestState();
        const vec2f& x = f.GetLocation();
        u_int s = u_int(f.GetScale());

        // if feature is still alive and we have a integral image at its scale but
        // it violates the distance assumption, kill it
        if((*it)->GetStatus() && s<imgs.size() && imgs[s].EvaluateApproximately<float>(x,hsize2)>1)
            (*it)->SetStatus(false);



    }

}

void CSimpleTracker::UpdateDescriptors(const std::vector<cv::Mat>& pyramid) {

    // smooth image for gradient computation
    /*vector<Mat> pyrsmooth;
    for(size_t s=0; s<pyramid.size(); s++) {

        Mat imsmooth;
        GaussianBlur(pyramid[s],imsmooth,Size(0,0),m_params->GetDoubleParameter("GRAD_SMOOTH_SIGMA"));

        pyrsmooth.push_back(imsmooth);

    }*/

    list<shared_ptr<mytracklet> >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus()) {

            imfeature& x = (*it)->GetLatestState();
            const vec2f& u0 = x.GetLocation();
            int s = int(x.GetScale());

            // create region of interest
            CRectangle<double> droi(u0.Get(0),
                                    u0.Get(1),
                                    m_params->GetIntParameter("DESCRIPTOR_HSIZE"),
                                    m_params->GetIntParameter("DESCRIPTOR_HSIZE"));

            // adjust region to scale
            if(s)
                droi.Scale(1.0/double(2<<(s-1)));

            // compute brief descriptor
            CBRIEF* briefdesc1 = new CBRIEF(droi);
            shared_ptr<CAbstractDescriptor> brief(briefdesc1);
            //brief->Compute(pyrsmooth[s]);
            brief->Compute(pyramid[s]);
            x.AttachDescriptor("BRIEF",brief);

            // interpret tracklet as simple tracker tracklet
            shared_ptr<CSimpleTrackerTracklet> tracklet = static_pointer_cast<CSimpleTrackerTracklet>(*it);
            float quality = 0;

            // check whether this has just been added or not
            if((*it)->GetCreationTime()==m_global_t)
                tracklet->m_reference_feature.AttachDescriptor("BRIEF",brief);
            else {

                if(tracklet->m_reference_feature.HasDescriptor("BRIEF")) {

                    shared_ptr<CAbstractDescriptor> desc0 = tracklet->m_reference_feature.GetDescriptor("BRIEF");
                    shared_ptr<CBRIEF> briefdesc0 = static_pointer_cast<CBRIEF>(desc0);
                    quality = briefdesc0->Distance(*briefdesc1);

                }

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
                //temp->Compute(pyrsmooth[s]);
                temp->Compute(pyramid[s]);

                shared_ptr<CAbstractDescriptor> idg(temp);
                x.AttachDescriptor("GRAD",idg);

            }

            if(m_params->GetIntParameter("COMPUTE_HOG")) {

                // compute HoG
                CHistogramOfGradients* temp = new CHistogramOfGradients(droi);

                //temp->Compute(imsmooth);
                temp->Compute(pyramid[s]);
                //temp->Normalize(2,1e-3);

                shared_ptr<CAbstractDescriptor> pdesc(temp);
                x.AttachDescriptor("HOG",pdesc);

            }

        }

    }

}

} // end of namespace
