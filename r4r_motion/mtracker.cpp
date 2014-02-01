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

CMotionTracker::CMotionTracker(CParameters* params, CPinholeCam &cam):
        CSimpleTracker(params),
        m_cam(cam),
        m_motion() {

    m_motion.push_back(vecf(6));

}

CView<float> CMotionTracker::GetLatestView() {

    vecf m0 = m_motion.back();

    CRigidMotion<float,3> F(m0);
    CView<float> view(m_cam,F);

    return view;

}

bool CMotionTracker::Update(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

    // do LK tracking at every time instance
    CSimpleTracker::Update(pyramid0,pyramid1);

    size_t kfr = (size_t)m_params->GetIntParameter("KEYFRAME_RATE");

    if(m_global_t>0 && (m_global_t)%kfr==0) {

        vecf m0 = m_motion.back();
        CRigidMotion<float,3> F0inv(m0);
        F0inv.Invert();

        vector<vec3f> xs;
        vector<vec2f> p0s, p1s, p1ss;
        vector<mytracklet*> trackletss2i, trackletsi2i;

        list<shared_ptr<mytracklet> >::iterator it;

        // collect image-to-image and scene-to-image correspondences
        for(it=m_data.begin(); it!=m_data.end(); it++) {

            vec2f p1 = (*it)->GetLatestLocationAtNativeScale();

            // cast to special tracklet type
            shared_ptr<CSimpleTrackerTracklet> tracklet = static_pointer_cast<CSimpleTrackerTracklet>(*it);

            if((*it)->GetStatus() && !tracklet->m_reference_feature.GetLocation().IsZero()) {

                // extract scene point
                xs.push_back(tracklet->m_reference_feature.GetLocation());
                p1ss.push_back(p1);
                trackletss2i.push_back((*it).get());

            }
            else if((*it)->GetStatus() && (*it)->GetLifetime()>(size_t)kfr) {

                vec2f u0 = (*it)->GetPastLocationAtNativeScale((size_t)kfr);

                p0s.push_back(u0);
                p1s.push_back(p1);
                trackletsi2i.push_back((*it).get());

            }

        }

        pair<vector<vec2f>,vector<vec2f> > corri2i(p0s,p1s);
        pair<vector<vec3f>,vector<vec2f> > corrs2i(xs,p1ss);

        cout << "No of correspondences (2d-2d/3d-2d): " << p0s.size() << " " << xs.size() << endl;

        // init linear solver
        CPreconditioner<smatf,float> M;
        CConjugateGradientMethodLeastSquares<smatf,float> solver(M,
                                                                 m_params->GetIntParameter("CGLS_NITER"),                                                                                     m_params->GetDoubleParameter("CGLS_EPS"),
                                                                 true);

        // init least-squares problem
        CMagicSfM problem(m_cam,corri2i,corrs2i,F0inv);

        // set up LM method
        CLevenbergMarquardt<smatf,float> lms(problem,solver,m_params->GetDoubleParameter("LM_LAMBDA"));

        // initialize
        vecf& model = problem.Get();

        for(size_t i=0; i<p0s.size(); i++)
            model(i) = m_params->GetDoubleParameter("INIT_DISTANCE");

        for(size_t i=0; i<6; i++)
            model(p0s.size()+i) = m0(i);

        // iterate
        vecf r = lms.Iterate(m_params->GetIntParameter("LM_NITER_OUTER"),
                             m_params->GetIntParameter("LM_NITER_INNER"),
                             m_params->GetDoubleParameter("LM_EPS"),
                             false,
                             false);

        // add new points to the map
        float threshold = m_params->GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_DEPTH");

        // create tentative map points
        vector<vec3f> x = m_cam.CAbstractCam::Normalize(corri2i.first);
        for(size_t i=0; i<x.size(); i++)
            x[i] = x[i]*model.Get(i);
        x = F0inv.Transform(x);

        // clear last map
        m_map.clear();

        for(size_t i=0; i<p0s.size(); i++) {

            if(fabs(r.Get(i))<threshold && fabs(r.Get(p0s.size()+i))<threshold) {

                // inject into map container
                m_map.push_back(pair<vec2f,vec3f>(corri2i.first[i],x[i]));

                // cast to special tracklet type and set reference point
                CSimpleTrackerTracklet* tracklet = reinterpret_cast<CSimpleTrackerTracklet*>(trackletsi2i[i]);
                tracklet->m_reference_feature.SetLocation(x[i]);

            }
            else
                trackletsi2i[i]->SetStatus(false);

        }

        // process motion
        threshold = m_params->GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_MOTION");

        vecf m1(6);
        for(size_t i=0; i<6; i++)
            m1(i) = model(p0s.size()+i);

        m_motion.push_back(m1);

        // reject outliers
        for(size_t i=0; i<xs.size(); i++) {

            if(fabs(r.Get(2*p0s.size()+i))>threshold || fabs(r.Get(2*p0s.size()+xs.size()+i))>threshold)
                trackletss2i[i]->SetStatus(false);

       }

    }

    return 0;

}

/*bool CMotionTracker::UpdateDescriptors(std::vector<cv::Mat>& pyramid) {

    list<shared_ptr<CCircularTracklet> >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus()) {

            imfeature& x = (*it)->GetLatestState();
            vec2f u0 = x.GetLocation();
            int s = int(x.GetScale());

            // create new feature
            CRectangle<double> droi(u0.Get(0),
                                    u0.Get(1),
                                    m_params->GetIntParameter("DESCRIPTOR_HSIZE"),
                                    m_params->GetIntParameter("DESCRIPTOR_HSIZE"));

            // adjust region to scale
            if(s)
                droi.Scale(1.0/double(2<<(s-1)));

            //droi.RotateTo(motion(5));

            // compute brief descriptor
            CBRIEF* briefdesc1 = new CBRIEF(droi);
            shared_ptr<CAbstractDescriptor> brief(briefdesc1);
            brief->Compute(pyramid[s]);
            //brief->Compute(pyramid[s]);
            x.AttachDescriptor("BRIEF",brief);

            // compute quality
            shared_ptr<CSimpleTrackerTracklet> tracklet = static_pointer_cast<CSimpleTrackerTracklet>(*it);
            float quality = 0;
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

        }

    }

    return 0;

}*/

CMagicSfM::CMagicSfM(CPinholeCam cam, pair<vector<vec2f>,vector<vec2f> >& corri2i, pair<vector<vec3f>,vector<vec2f> >& corrs2i, CRigidMotion<float,3> F0inv):
    CLeastSquaresProblem<smatf,float>::CLeastSquaresProblem(2*(corri2i.first.size()+corrs2i.first.size()),corri2i.first.size()+6),
	m_cam(cam),
	m_corri2i(corri2i),
	m_corrs2i(corrs2i),
	m_F0inv(F0inv) {

}

void CMagicSfM::ComputeResidual(vecf& r) {

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
    for(size_t i=0; i<m; i++) {

        // projection error
        vec2f dp = p1p[i] - m_corri2i.second[i];

        // copy weighted error
        r(i)   = m_weights.Get(i)*dp.Get(0);
        r(m+i) = m_weights.Get(m+i)*dp.Get(1);

    }

    // project map points to second view, they are in world coordinates
    p1p = view1.Project(m_corrs2i.first);

    // copy errors
    for(size_t i=0; i<n; i++) {

        vec2f dp = p1p[i] - m_corrs2i.second[i];

        // set residual
        r(2*m+i)   = m_weights.Get(2*m+i)*dp.Get(0);
        r(2*m+n+i) = m_weights.Get(2*m+n+i)*dp.Get(1);

    }

}

void CMagicSfM::ComputeResidualAndJacobian(vecf& r, smatf& J) {

    /* Jacobian w.r.t.
     * - rotation need backprojected point in world coordinates,
     * - depths need viewing directions of frame 0 in frame 1 coordinates,
     * - translations only need Jacobian of pinhole projection into second view,
     * - all involve Jacobian of second pinhole projection.
     *
    */

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


    // create second view
    CView<float> view1(m_cam,F1);

    // normalize pixel in first frame
    vector<vec3f> x0n = m_cam.CAbstractCam::Normalize(m_corri2i.first);

    // multiply by current depth
    vector<vec3f> x0;
    for(size_t i=0; i<x0n.size(); i++)
        x0.push_back(x0n[i]*m_model.Get(i));

    // transform point to world coordinates, keep this for Jacobian w.r.t. to rotation
    x0 = m_F0inv.Transform(x0);

    // transform viewing direction from one frame to the other (depth derivative)
    x0n = Fr.DifferentialTransform(x0n);

    // projection error for image-to-image correspondences
    for(size_t i=0; i<m; i++) {

        float wi, wim;
        wi = m_weights.Get(i);
        wim = m_weights.Get(m+i);

        // transform point into frame 1
        vec3f x1 = F1.Transform(x0[i]);

        // project and compute Jacobian
        vec2f p1p;
        matf Jpi(2,3);
        m_cam.Project(x1,p1p,Jpi);

        // projection error
        vec2f dp = p1p - m_corri2i.second[i];

        // set residual
        r(i)   = wi*dp.Get(0);
        r(m+i) = wim*dp.Get(1);

        // depth derivatives
        J(i,i)   = wi*(Jpi.Get(0,0)*x0n[i].Get(0) + Jpi.Get(0,2)*x0n[i].Get(2));
        J(m+i,i) = wim*(Jpi.Get(1,1)*x0n[i].Get(1) + Jpi.Get(1,2)*x0n[i].Get(2));

        // translational derivative
        J(i,m)   = wi*Jpi(0,0);     J(i,m+1)   = wi*Jpi(0,1);	J(i,m+2)   = wi*Jpi(0,2);
        J(m+i,m) = wim*Jpi(1,0); 	J(m+i,m+1) = wim*Jpi(1,1); 	J(m+i,m+2) = wim*Jpi(1,2);

        // rotational derivatives
        vec3f do1 = DR1x*x0[i];
        vec3f do2 = DR1y*x0[i];
        vec3f do3 = DR1z*x0[i];

        J(i,m+3)   = wi*(Jpi.Get(0,0)*do1.Get(0) + Jpi.Get(0,2)*do1.Get(2));
        J(m+i,m+3) = wim*(Jpi.Get(1,1)*do1.Get(1) + Jpi.Get(1,2)*do1.Get(2));

        J(i,m+4)   = wi*(Jpi.Get(0,0)*do2.Get(0) + Jpi.Get(0,2)*do2.Get(2));
        J(m+i,m+4) = wim*(Jpi.Get(1,1)*do2.Get(1) + Jpi.Get(1,2)*do2.Get(2));

        J(i,m+5)   = wi*(Jpi.Get(0,0)*do3.Get(0) + Jpi.Get(0,2)*do3.Get(2));
        J(m+i,m+5) = wim*(Jpi.Get(1,1)*do3.Get(1) + Jpi.Get(1,2)*do3.Get(2));

    }

    // scene-to-image correspondences
    for(size_t i=0; i<n; i++) {

        float w2mi, w2min;
        w2mi = m_weights.Get(2*m+i);
        w2min = m_weights.Get(2*m+n+i);

        // transform map points to frame 1
        vec3f x1 = F1.Transform(m_corrs2i.first[i]);

        // projection into second image
        vec2f p1p;
        matf Jpi(2,3);
        m_cam.Project(x1,p1p,Jpi);

        // projection error
        vec2f dp = p1p - m_corrs2i.second[i];

        // set residual
        r(2*m+i)   = w2mi*dp.Get(0);
        r(2*m+n+i) = w2min*dp.Get(1);

        // translational derivative
        J(2*m+i,m)   = w2mi*Jpi(0,0);   J(2*m+i,m+1)   = w2mi*Jpi(0,1);		J(2*m+i,m+2)   = w2mi*Jpi(0,2);
        J(2*m+n+i,m) = w2min*Jpi(1,0); 	J(2*m+n+i,m+1) = w2min*Jpi(1,1);    J(2*m+n+i,m+2) = w2min*Jpi(1,2);

        // rotational derivatives
        vec3f do1 = DR1x*m_corrs2i.first[i];
        vec3f do2 = DR1y*m_corrs2i.first[i];
        vec3f do3 = DR1z*m_corrs2i.first[i];

        J(2*m+i,m+3)   = w2mi*(Jpi.Get(0,0)*do1.Get(0) + Jpi.Get(0,2)*do1.Get(2));
        J(2*m+n+i,m+3) = w2min*(Jpi.Get(1,1)*do1.Get(1) + Jpi.Get(1,2)*do1.Get(2));

        J(2*m+i,m+4)   = w2mi*(Jpi.Get(0,0)*do2.Get(0) + Jpi.Get(0,2)*do2.Get(2));
        J(2*m+n+i,m+4) = w2min*(Jpi.Get(1,1)*do2.Get(1) + Jpi.Get(1,2)*do2.Get(2));

        J(2*m+i,m+5)   = w2mi*(Jpi.Get(0,0)*do3.Get(0) + Jpi.Get(0,2)*do3.Get(2));
        J(2*m+n+i,m+5) = w2min*(Jpi.Get(1,1)*do3.Get(1) + Jpi.Get(1,2)*do3.Get(2));

    }

}

vecf CMagicSfM::ComputeDispersion(const vecf& r) {

    vecf sigma(r.NElems());

    size_t m = m_corri2i.first.size();
    size_t n = m_corrs2i.first.size();

    // get pointer to the data
    float* rdata = r.Data().get();

    // copy into subarrays (one for each dimension and part of the functional)
    vecf ri2iu(m,shared_ptr<float>(&rdata[0]));
    vecf ri2iv(m,shared_ptr<float>(&rdata[m]));
    vecf rs2iu(n,shared_ptr<float>(&rdata[2*m]));
    vecf rs2iv(n,shared_ptr<float>(&rdata[2*m+n]));

    float si2iu = 1.0/(1.4826*ri2iu.MAD());
    float si2iv = 1.0/(1.4826*ri2iv.MAD());
    float ss2iu = 1.0/(1.4826*rs2iu.MAD());
    float ss2iv = 1.0/(1.4826*rs2iv.MAD());

    for(size_t i=0; i<m; i++) {

        sigma(i) = si2iu;
        sigma(m+i) = si2iv;

    }

    for(size_t i=0; i<n; i++) {

        sigma(2*m+i) = ss2iu;
        sigma(2*m+n+i) = ss2iv;

    }

    return sigma;

}

} // end of namespace
