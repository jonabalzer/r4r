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

#include "mtracker.h"
#include "lm.h"
#include "descriptor.h"
#include "rutils.h"
#include "trafo.h"
#include "interp.h"

#include <fstream>

using namespace cv;

namespace R4R {

CMotionTracker::CMotionTracker(const CParameters* params, const CPinholeCam<float>& cam):
        CSimpleTracker(params),
        m_cam(cam),
        m_motion(),
        m_map(1000) {

    // set initial frame to the identity
    m_motion.Update(0,0,0,0,0,0);

}

CView<float> CMotionTracker::GetLatestView() {

    CRigidMotion<float,3> F = m_motion.GetLatestState();
    return CView<float>(m_cam,F);

}

void CMotionTracker::Update(const vector<Mat>& pyramid0, const vector<Mat>& pyramid1) {

    // do LK tracking at every time instance
    CSimpleTracker::Update(pyramid0,pyramid1);

    size_t kfr = (size_t)m_params->GetIntParameter("KEYFRAME_RATE");

    if(m_global_t>0 && m_global_t%kfr==0) {

        CRigidMotion<float,3> F0inv = m_motion.GetLatestState();
        F0inv.Invert();

        vector<vec3f> xs;                                   // map points of view 0
        vector<vec2f> p0s, p1s, p1ss;                       //
        vector<mytracklet*> trackletss2i, trackletsi2i;

        list<shared_ptr<mytracklet> >::iterator it;

        // collect image-to-image and scene-to-image correspondences
        for(it=m_data.begin(); it!=m_data.end(); it++) {

            vec2f p1 = (*it)->GetLatestLocationAtNativeScale();

            // cast to special tracklet type
            shared_ptr<CMotionTrackerTracklet> tracklet = static_pointer_cast<CMotionTrackerTracklet>(*it);

            if((*it)->GetStatus() && tracklet->m_has_point) {

                // extract scene point
                xs.push_back(tracklet->m_pmap_point->GetLocation());
                p1ss.push_back(p1);
                trackletss2i.push_back((*it).get());

            } // only add more im2im correspondences at key frames
            else if((*it)->GetStatus() && (*it)->GetLifetime()>(size_t)kfr) {

                // check for the ringbuffer wether this is ok!!!
                vec2f u0 = (*it)->GetPastLocationAtNativeScale((size_t)kfr);

                p0s.push_back(u0);
                p1s.push_back(p1);
                trackletsi2i.push_back((*it).get());

            }

        }

        cout << "No of correspondences (2d-2d/3d-2d): " << p0s.size() << " " << xs.size() << endl;

        // nothing to do?
        if(p0s.size()==0 && xs.size()==0)
            return;

        // group into pair
        pair<vector<vec2f>,vector<vec2f> > corri2i(p0s,p1s);
        pair<vector<vec3f>,vector<vec2f> > corrs2i(xs,p1ss);

        // init linear solver
        CPreconditioner<CCSRMatrix<float>,float> M;
        CConjugateGradientMethodLeastSquares<CCSRMatrix<float>,float> solver(M,
                                                                             m_params->GetIntParameter("CGLS_NITER"),                                                                                     m_params->GetDoubleParameter("CGLS_EPS"),
                                                                             true);

        // init least-squares problem
        CMagicSfM problem(m_cam,corri2i,corrs2i,F0inv);

        // set up LM method
        CLevenbergMarquardt<CCSRMatrix<float>,float> lms(problem,solver,m_params->GetDoubleParameter("LM_LAMBDA"));

        // initialize
        vecf& model = problem.Get();

        for(size_t i=0; i<p0s.size(); i++)
            model(i) = m_params->GetDoubleParameter("INIT_DISTANCE");

        // get motion
        const list<CVector<float,6> >& mdata = m_motion.GetData();
        CVector<float,6> m0 = mdata.back();
        for(size_t i=0; i<6; i++)
            model(p0s.size()+i) = m0.Get(i);

        // iterate
        CBiSquareWeightFunction<float> weight;
        vecf r = lms.Iterate(m_params->GetIntParameter("LM_NITER_OUTER"),
                             weight,
                             m_params->GetIntParameter("LM_NITER_INNER"),
                             m_params->GetDoubleParameter("LM_EPS"),
                             true,
                             true);

        // add new points to the map
        float threshold = m_params->GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_DEPTH");

        // create tentative map points
        vector<vec3f> x = m_cam.CAbstractCam::Normalize(corri2i.first);
        for(size_t i=0; i<x.size(); i++)
            x[i] = x[i]*model.Get(i);
        x = F0inv.Transform(x);

        size_t row = 0;
        for(size_t i=0; i<p0s.size(); i++) {

            if(fabs(r.Get(row))<threshold && fabs(r.Get(row+1))<threshold) {

                // cast to special tracklet type and set reference point
                CMotionTrackerTracklet* tracklet = reinterpret_cast<CMotionTrackerTracklet*>(trackletsi2i[i]);

                // inject into map container
                CInterestPoint<float,3> feature(x[i],
                                                tracklet->GetLatestState().GetScale(),
                                                tracklet->GetLatestState().GetQuality());
                tracklet->m_pmap_point = m_map.Insert(feature);
                tracklet->m_has_point = true;

            }
            else
                trackletsi2i[i]->SetStatus(false);

            row += 2;

        }

        // process motion
        threshold = m_params->GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_MOTION");

        // inject into trajectory
        m_motion.Update(model.Get(p0s.size()),
                        model.Get(p0s.size()+1),
                        model.Get(p0s.size()+2),
                        model.Get(p0s.size()+3),
                        model.Get(p0s.size()+4),
                        model.Get(p0s.size()+5));

        // reject outliers
        for(size_t i=0; i<xs.size(); i++) {

            if(fabs(r.Get(row))>threshold || fabs(r.Get(row+1))>threshold)
                trackletss2i[i]->SetStatus(false);

            row += 2;

       }

    }

}

void CMotionTracker::AddTracklets(const std::vector<Mat>& pyramid) {

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

                            // create new tracklet with feature, set the iterator to the end of the map
                            CMotionTrackerTracklet* tracklet = new CMotionTrackerTracklet(m_global_t,x,m_params->GetIntParameter("BUFFER_LENGTH"));
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

CMagicSfM::CMagicSfM(CPinholeCam<float> cam, pair<vector<vec2f>,vector<vec2f> >& corri2i, pair<vector<vec3f>,vector<vec2f> >& corrs2i, CRigidMotion<float,3> F0inv):
    CLeastSquaresProblem<CCSRMatrix<float>,float>::CLeastSquaresProblem(2*(corri2i.first.size()+corrs2i.first.size()),corri2i.first.size()+6),
	m_cam(cam),
	m_corri2i(corri2i),
	m_corrs2i(corrs2i),
	m_F0inv(F0inv) {

}

void CMagicSfM::ComputeResidual(vecf& r) const {

    size_t m = m_corri2i.first.size();
    size_t n = m_corrs2i.first.size();

    // actual transformation from world to second frame
    CRigidMotion<float,3> F1(m_model.Get(m),
                             m_model.Get(m+1),
                             m_model.Get(m+2),
                             m_model.Get(m+3),
                             m_model.Get(m+4),
                             m_model.Get(m+5));

    // create second view
    CView<float> view1(m_cam,F1);

    // normalize pixel in first frame, keep this in Jacobian for flow
    vector<vec3f> x = m_cam.CAbstractCam::Normalize(m_corri2i.first);

    // multiply by current depth
    for(size_t i=0; i<x.size(); i++)
        x[i] = x[i]*m_model.Get(i);

    // transform to world coordinates
    x = m_F0inv.Transform(x);

    // second view
    vector<vec2f> p1p = view1.Project(x);

    // enter error in residual vector
    size_t row = 0;
    for(size_t i=0; i<m; i++) {

        // projection error
        vec2f dp = p1p[i] - m_corri2i.second[i];

        // copy weighted error
        r(row) = dp.Get(0);
        row++;
        r(row) = dp.Get(1);
        row++;

    }

    // project map points to second view, they are in world coordinates
    p1p = view1.Project(m_corrs2i.first);

    // copy errors
    for(size_t i=0; i<n; i++) {

        vec2f dp = p1p[i] - m_corrs2i.second[i];

        // set residual
        r(row) = dp.Get(0);
        row++;
        r(row) = dp.Get(1);
        row++;

    }

}

void CMagicSfM::ComputeResidualAndJacobian(vecf& r, CCSRMatrix<float> &J) const {

    /* Jacobian w.r.t.
     * - rotation need backprojected point in world coordinates,
     * - depths need viewing directions of frame 0 in frame 1 coordinates,
     * - translations only need Jacobian of pinhole projection into second view,
     * - all involve Jacobian of second pinhole projection.
     *
    */

    // init matrix data, the sizes are set by the calling LM routine
    vector<size_t>* rowptr(new vector<size_t>());
    vector<size_t>* cols(new vector<size_t>());
    vector<float>* vals(new vector<float>());

    // nnz counter
    size_t nnz = 0;
    rowptr->push_back(nnz);

    // get sizes for easier book-keeping of residual indices
    size_t m = m_corri2i.first.size();
    size_t n = m_corrs2i.first.size();

    // actual transformation from world to second frame
    CRigidMotion<float,3> F1(m_model.Get(m),
                             m_model.Get(m+1),
                             m_model.Get(m+2),
                             m_model.Get(m+3),
                             m_model.Get(m+4),
                             m_model.Get(m+5));


    // concatenation of F1*m_F0inv, only Jacobian of this is needed later
    CTransformation<float,3> F = F1*m_F0inv;
    CRigidMotion<float,3> Fr = reinterpret_cast<CRigidMotion<float,3>& >(F);

    // derivatives of *rotations* in DF1
    matf DR1x(3,3), DR1y(3,3), DR1z(3,3);
    CDifferentialRotation<float,3>::Rodrigues(m_model.Get(m+3),
                                              m_model.Get(m+4),
                                              m_model.Get(m+5),
                                              nullptr,
                                              DR1x.Data().get(),
                                              DR1y.Data().get(),
                                              DR1z.Data().get());

    // normalize pixel in first frame
    vector<vec3f> x0n = m_cam.CAbstractCam::Normalize(m_corri2i.first);

    // multiply by current depth
    vector<vec3f> x0;
    x0.reserve(x0n.size());
    for(size_t i=0; i<x0n.size(); i++)
        x0.push_back(x0n[i]*m_model.Get(i));

    // transform point to world coordinates, keep this for Jacobian w.r.t. to rotation
    x0 = m_F0inv.Transform(x0);

    // transform viewing direction from one frame to the other (depth derivative)
    x0n = Fr.DifferentialTransform(x0n);

    // projection error for image-to-image correspondences
    size_t row = 0;
    size_t mp1, mp2, mp3, mp4, mp5;
    mp1 = m + 1;
    mp2 = m + 2;
    mp3 = m + 3;
    mp4 = m + 4;
    mp5 = m + 5;

    for(size_t i=0; i<m; i++) {

        // transform point into frame 1
        vec3f x1 = F1.Transform(x0[i]);

        // project and compute Jacobian
        vec2f p1p;
        matf Jpi(2,3);
        m_cam.Project(x1,p1p,Jpi);

        // projection error
        vec2f dp = p1p - m_corri2i.second[i];

        // rotational derivatives
        vec3f do1 = DR1x*x0[i];
        vec3f do2 = DR1y*x0[i];
        vec3f do3 = DR1z*x0[i];

        // one row for the u coordinate
        float wi = m_weights.Get(row);
        r(row) = wi*dp.Get(0);

        // depth derivatives
        vals->push_back(wi*(Jpi.Get(0,0)*x0n[i].Get(0) + Jpi.Get(0,2)*x0n[i].Get(2)));
        cols->push_back(i);

        // translational derivative
        vals->push_back(wi*Jpi.Get(0,0));
        vals->push_back(wi*Jpi.Get(0,1));
        vals->push_back(wi*Jpi.Get(0,2));
        cols->push_back(m);
        cols->push_back(mp1);
        cols->push_back(mp2);

        // rotational derivative
        vals->push_back(wi*(Jpi.Get(0,0)*do1.Get(0) + Jpi.Get(0,2)*do1.Get(2)));
        vals->push_back(wi*(Jpi.Get(0,0)*do2.Get(0) + Jpi.Get(0,2)*do2.Get(2)));
        vals->push_back(wi*(Jpi.Get(0,0)*do3.Get(0) + Jpi.Get(0,2)*do3.Get(2)));
        cols->push_back(mp3);
        cols->push_back(mp4);
        cols->push_back(mp5);

        // increment row and nnz counter
        row++;
        nnz += 7;
        rowptr->push_back(nnz);

        // same procedure for the v coordinate
        wi = m_weights.Get(row);
        r(row) = wi*dp.Get(1);

        // depth derivative
        vals->push_back(wi*(Jpi.Get(1,1)*x0n[i].Get(1) + Jpi.Get(1,2)*x0n[i].Get(2)));
        cols->push_back(i);

        // translational derivative
        vals->push_back(wi*Jpi(1,0));
        vals->push_back(wi*Jpi(1,1));
        vals->push_back(wi*Jpi(1,2));
        cols->push_back(m);
        cols->push_back(mp1);
        cols->push_back(mp2);

        // rotational derivative
        vals->push_back(wi*(Jpi.Get(1,1)*do1.Get(1) + Jpi.Get(1,2)*do1.Get(2)));
        vals->push_back(wi*(Jpi.Get(1,1)*do2.Get(1) + Jpi.Get(1,2)*do2.Get(2)));
        vals->push_back(wi*(Jpi.Get(1,1)*do3.Get(1) + Jpi.Get(1,2)*do3.Get(2)));
        cols->push_back(mp3);
        cols->push_back(mp4);
        cols->push_back(mp5);

        // update counters
        row++;
        nnz += 7;
        rowptr->push_back(nnz);

    }

    // scene-to-image correspondences
    for(size_t i=0; i<n; i++) {

        // transform map points to frame 1
        vec3f x1 = F1.Transform(m_corrs2i.first[i]);

        // projection into second image
        vec2f p1p;
        matf Jpi(2,3);
        m_cam.Project(x1,p1p,Jpi);

        // projection error
        vec2f dp = p1p - m_corrs2i.second[i];

        // rotational derivatives
        vec3f do1 = DR1x*m_corrs2i.first[i];
        vec3f do2 = DR1y*m_corrs2i.first[i];
        vec3f do3 = DR1z*m_corrs2i.first[i];

        // first the u direction
        float wi = m_weights.Get(row);
        r(row) = wi*dp.Get(0);

        // translational derivative
        vals->push_back(wi*Jpi.Get(0,0));
        vals->push_back(wi*Jpi.Get(0,1));
        vals->push_back(wi*Jpi.Get(0,2));
        cols->push_back(m);
        cols->push_back(mp1);
        cols->push_back(mp2);

        // rotational derivatives
        vals->push_back(wi*(Jpi.Get(0,0)*do1.Get(0) + Jpi.Get(0,2)*do1.Get(2)));
        vals->push_back(wi*(Jpi.Get(0,0)*do2.Get(0) + Jpi.Get(0,2)*do2.Get(2)));
        vals->push_back(wi*(Jpi.Get(0,0)*do3.Get(0) + Jpi.Get(0,2)*do3.Get(2)));
        cols->push_back(mp3);
        cols->push_back(mp4);
        cols->push_back(mp5);

        // increment row and nnz counter
        row++;
        nnz += 6;
        rowptr->push_back(nnz);

        // now the v direction
        wi = m_weights.Get(row);
        r(row) = wi*dp.Get(1);

        // translational derivative
        vals->push_back(wi*Jpi.Get(1,0));
        vals->push_back(wi*Jpi.Get(1,1));
        vals->push_back(wi*Jpi.Get(1,2));
        cols->push_back(m);
        cols->push_back(mp1);
        cols->push_back(mp2);

        // rotational derivative
        vals->push_back(wi*(Jpi.Get(1,1)*do1.Get(1) + Jpi.Get(1,2)*do1.Get(2)));
        vals->push_back(wi*(Jpi.Get(1,1)*do2.Get(1) + Jpi.Get(1,2)*do2.Get(2)));
        vals->push_back(wi*(Jpi.Get(1,1)*do3.Get(1) + Jpi.Get(1,2)*do3.Get(2)));
        cols->push_back(mp3);
        cols->push_back(mp4);
        cols->push_back(mp5);

        // increment row and nnz counter
        row++;
        nnz += 6;
        rowptr->push_back(nnz);

    }

    J.SetData(rowptr,cols,vals);

}

} // end of namespace
