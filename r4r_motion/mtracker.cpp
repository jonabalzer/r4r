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
#include "basic.h"
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
        vector<CTracklet*> trackletss2i, trackletsi2i;

        list<shared_ptr<CTracklet> >::iterator it;

        // collect image-to-image and scene-to-image correspondences
        for(it=begin(); it!=end(); it++) {

            vec2f p1 = (*it)->GetLatestLocationAtNativeScale();

            if((*it)->GetStatus() && (*it)->front().HasDescriptor("3DPOINT")) {

                // get access to the scene point
                shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("3DPOINT");
                CDescriptor<vecf>* cdesc = (CDescriptor<vecf>*)desc.get();
                vec3f x(cdesc->Get());

                xs.push_back(x);
                p1ss.push_back(p1);
                trackletss2i.push_back((*it).get());

            }
            else if((*it)->GetStatus() && (*it)->size()>(size_t)kfr) {

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
        CConjugateGradientMethodLeastSquares<smatf,float> solver(CPreconditioner<smatf,float>(),
                                                                 m_params->GetIntParameter("CGLS_NITER"),                                                                                     m_params->GetDoubleParameter("CGLS_EPS"),                                                                                     true);

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
                             true);

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

                // get pointer to initial feature
                shared_ptr<CAbstractDescriptor> descx(new CDescriptor<vecf>(vecf(x[i])));
                trackletsi2i[i]->front().AttachDescriptor("3DPOINT",descx);

                // also store initial depths for descriptor canonization
//                vec z(1);
//                z(0) = x(i);
//                shared_ptr<CDescriptor> descz(new CDescriptor<vec>(z));
//                feat->AttachDescriptor("INITDEPTH",descz);

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

bool CMotionTracker::UpdateDescriptors(std::vector<cv::Mat>& pyramid) {

    // get actual pose
    /*vec motion = m_motion.back();
    CRigidMotion<float,3> F(motion);*/

    list<shared_ptr<CTracklet> >::iterator it;

    for(it=begin(); it!=end(); it++) {

        if((*it)->GetStatus()) {

            imfeature& x = (*it)->GetLatestState();
            vec2f u0 = x.GetLocation();
            u_int s = u_int(x.GetScale()+0.5);

            // create new feature
            CRectangle<double> droi(u0.Get(0),
                                    u0.Get(1),
                                    m_params->GetIntParameter("DESCRIPTOR_HSIZE"),
                                    m_params->GetIntParameter("DESCRIPTOR_HSIZE"));

            // adjust region to scale
            droi.Scale(1.0/pow(2,s));
            //droi.RotateTo(motion(5));

            /*// see the feature has an initial depth
            if((*it)->front()->HasDescriptor("INITDEPTH")) {

                // get initial depth
                shared_ptr<CDescriptor> descz = (*it)->front()->GetDescriptor("INITDEPTH");
                CDescriptor<vec>* cdescz = (CDescriptor<vec>*)descz.get();
                vec z = cdescz->Get();

                // get point in local coordinates
                shared_ptr<CDescriptor> descx = (*it)->front()->GetDescriptor("3DPOINT");
                CDescriptor<vec>* cdescx = (CDescriptor<vec>*)descx.get();
                vec pt = CLinearAlgebra::TransformPoint(F, cdescx->Get());

                // scale by z0/z
                droi.Scale(z(0)/pt(2));

            }*/

            // compute brief descriptor
            CBRIEF* briefdesc1 = new CBRIEF(droi);
            shared_ptr<CAbstractDescriptor> brief(briefdesc1);
            brief->Compute(pyramid[s]);
            x.AttachDescriptor("BRIEF",brief);

            // compute quality
            shared_ptr<CAbstractDescriptor> desc0 = (*it)->front().GetDescriptor("BRIEF");
            double quality = 1;

            if(desc0!=nullptr) {

                shared_ptr<CBRIEF> briefdesc0 = static_pointer_cast<CBRIEF>(desc0);
                quality = briefdesc0->Distance(*briefdesc1);

            }

            x.SetQuality(quality);

        }

    }

    return 0;

}

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

    // projection error for image-to-image correspondences
/*	for(size_t i=0; i<m_corri2i.size(); i++) {

		// pixel in img 0 and 1
        vec2f p0, p1;
        p0 = m_corri2i[i].first;
        p1 = m_corri2i[i].second;

		// estimate of the scene point x0 in camera system of first frame
        vec3f x0h = m_cam.Normalize(p0);
        vec3f x0 = x0h*m_model.Get(i);
        x0 = m_F0inv.Transform(x0);

		// point transformed from world to local camera system of second frame
        vec3f x1 = F1.Transform(x0);

		// projection of x0 into second image
        vec2f p1p = m_cam.Project(x1);

		// projection error
        vec2f dp = p1p - p1;

        // set residual
        r(i) = m_weights.Get(i)*dp.Get(0);
        r(m_corri2i.size()+i) = m_weights.Get(m_corri2i.size()+i)*dp.Get(1);

	}

	// projection error for scene-to-image correspondences
	for(size_t i=0; i<m_corrs2i.size(); i++) {

        vec3f x0 = m_corrs2i[i].first;		// scene point in world coordinates
        vec2f p1 = m_corrs2i[i].second;     // projection into second image

		// rigid transform from wc to second cam
        vec3f x1 = F1.Transform(x0);

        // predicted projection of x1 into second image
        vec2f p1p = m_cam.Project(x1);

		// projection error
        vec2f dp = p1p - p1;

		// set residual
        r(2*m_corri2i.size()+i) = m_weights.Get(2*m_corri2i.size()+i)*dp.Get(0);
        r(2*m_corri2i.size()+m_corrs2i.size()+i) = m_weights.Get(2*m_corri2i.size()+m_corrs2i.size()+i)*dp.Get(1);

    }*/

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

//CSfMTrackerUpdate::CSfMTrackerUpdate(CCam cam, vector<CTracklet*> map, cv::Mat& img0, cv::Mat& img1, size_t hsize):
//	CLeastSquaresProblem(map.size()*(2*hsize+1)*(2*hsize+1),6),
//	m_cam(cam),
//	m_map(map),
//	m_img0(img0),
//	m_img1(img1),
//	m_hsize(hsize) {

//}

//void CSfMTrackerUpdate::ComputeResidualAndJacobian(vec& r, mat& J) {

//	// some sizes
//	size_t n = m_map.size();					// no of anchor pixels
//	size_t w = 2*m_hsize+1; 					// width of neighborhood
//	size_t m = w*w;								// no of pixels around each anchor

//	// actual transformation
//	vec t(3), omega(3);
//	for(size_t i=0; i<3; i++) {

//		t(i) = m_model(i);
//		omega(i) = m_model(3+i);

//	}

//	mat R(3,3), DRx(3,3), DRy(3,3), DRz(3,3);
//	CLinearAlgebra::Rodrigues(omega,R,DRx,DRy,DRz);

//	// projection errors and their derivatives
//	for(size_t i=0; i<n; i++) {

//		// anchor point in first frame
//		vec u0 = m_map[i]->GetLatestLocation();

//		// world point
//        shared_ptr<CAbstractDescriptor> desc = m_map[i]->front().GetDescriptor("3DPOINT");
//        CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
//		vec x0 = cdesc->Get();

//		// rigid transform from wc to second cam
//		vec x1 = R*x0 + t;

//		// projection of x1 into second image
//		vec u1 = m_cam.ProjectLocal(x1);

//		// init Jacobian w.r.t to local coordinates of projection into second frame
//		mat Jpi(2,3);
//		Jpi(0,0) = m_cam.m_f[0]/x1(2);
//		Jpi(1,1) = m_cam.m_f[1]/x1(2);

//		// every pixel in neighborhood of anchor contributes to energy
//		for(int j=-(int)m_hsize; j<=(int)m_hsize; j++) {

//			for(int k=-(int)m_hsize; k<=(int)m_hsize; k++) {

//				// compute row index
//				int row = m*i+w*(j+(int)m_hsize)+(int)m_hsize+k;

//				// interpolated intensities and gradients
//				double I0, I1, I1u, I1v;
//				I0 = CImageInterpolation::Bilinear(m_img0,u0(0)+j,u0(1)+k);
//				I1 = CImageInterpolation::Bilinear(m_img1,u1(0)+j,u1(1)+k);
//				I1u = CImageInterpolation::Gradient(m_img1,u1(0)+j,u1(1)+k,true);
//				I1v = CImageInterpolation::Gradient(m_img1,u1(0)+j,u1(1)+k,false);

//				// residual
//				r(row) = I1 - I0;

//				// update projection Jacobian
//				Jpi(0,2) = -(u1(0)+j-m_cam.m_c[0])/x1(2);
//				Jpi(1,2) = -(u1(1)+k-m_cam.m_c[1])/x1(2);

//				// translational part
//				J(row,0) = I1u*Jpi(0,0) + I1v*Jpi(1,0);
//				J(row,1) = I1u*Jpi(0,1) + I1v*Jpi(1,1);
//				J(row,2) = I1u*Jpi(0,2) + I1v*Jpi(1,2);

//				// rotational part
//				vec do1 = Jpi*(DRx*x0);
//				vec do2 = Jpi*(DRy*x0);
//				vec do3 = Jpi*(DRz*x0);
//				J(row,3) = I1u*do1(0) + I1v*do1(1);
//				J(row,4) = I1u*do2(0) + I1v*do2(1);
//				J(row,5) = I1u*do3(0) + I1v*do3(1);

//			}

//		}

//	}

//}

//void CSfMTrackerUpdate::ComputeResidual(vec& r) {
//	// some sizes
//	size_t n = m_map.size();					// no of anchor pixels
//	size_t w = 2*m_hsize+1; 					// width of neighborhood
//	size_t m = w*w;								// no of pixels around each anchor

//	// actual transformation
//	vec t(3), omega(3);
//	for(size_t i=0; i<3; i++) {

//		t(i) = m_model(i);
//		omega(i) = m_model(3+i);

//	}

//	mat R(3,3), DRx(3,3), DRy(3,3), DRz(3,3);
//	CLinearAlgebra::Rodrigues(omega,R,DRx,DRy,DRz);

//	// projection errors and their derivatives
//#pragma omp parallel for
//	for(size_t i=0; i<n; i++) {

//		// anchor point in first frame
//		vec u0 = m_map[i]->GetLatestLocation();

//		// world point
//        shared_ptr<CAbstractDescriptor> desc = m_map[i]->front().GetDescriptor("3DPOINT");
//        CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
//		vec x0 = cdesc->Get();

//		// rigid transform from wc to second cam
//		vec x1 = R*x0 + t;

//		// projection of x1 into second image
//		vec u1 = m_cam.ProjectLocal(x1);

//		// every pixel in neighborhood of anchor contributes to energy
//		for(int j=-(int)m_hsize; j<=(int)m_hsize; j++) {

//			for(int k=-(int)m_hsize; k<=(int)m_hsize; k++) {

//				// compute row index
//				int row = m*i+w*(j+(int)m_hsize)+(int)m_hsize+k;

//				// interpolated intensities and gradients
//				double I0, I1;
//				I0 = CImageInterpolation::Bilinear(m_img0,u0(0)+j,u0(1)+k);	// FIXME: do this outside
//				I1 = CImageInterpolation::Bilinear(m_img1,u1(0)+j,u1(1)+k);

//				// residual
//				r(row) = I1 - I0;


//			}

//		}

//	}

//}


}
