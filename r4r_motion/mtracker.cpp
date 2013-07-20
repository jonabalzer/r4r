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

map<CTracklet*,vec> CMotionTracker::GetMap() {

    map<CTracklet*,vec> pts;
    list<shared_ptr<CTracklet> >::iterator it;

    for(it=begin(); it!=end(); it++) {

        vec ptcolor(6);

        // check whether we have 3d data for feature in first frame
        if((*it)->front().HasDescriptor("3DPOINT")) {

            // get 3d point
            shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("3DPOINT");
            CDescriptor<vecf>* cdesc = (CDescriptor<vecf>*)desc.get();
            vecf pt = cdesc->Get();

            // copy data
            ptcolor(0) = pt.Get(0);
            ptcolor(1) = pt.Get(1);
            ptcolor(2) = pt.Get(2);

        }

        // check whether the point is colored
        if((*it)->front().HasDescriptor("COLOR")) {

            // color
            shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("COLOR");
            CDescriptor<vecf>* cdesc = (CDescriptor<vecf>*)desc.get();
            vecf color = cdesc->Get();

            // copy data
            ptcolor(3) = color.Get(0);
            ptcolor(4) = color.Get(1);
            ptcolor(5) = color.Get(2);

        } else {

            ptcolor(3) = 255;
            ptcolor(4) = 255;
            ptcolor(5) = 255;

        }

        // only attach if it has at least a point
        if((*it)->front().HasDescriptor("3DPOINT"))
            pts.insert(pair<CTracklet*,vec>((*it).get(),ptcolor));

    }

    return pts;

}

bool CMotionTracker::Update(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

    // do LK tracking at every time instance, FIXME: how to do descriptor aggregation with 3d info
    CSimpleTracker::Update(pyramid0,pyramid1);

    size_t kfr = (size_t)m_params->GetIntParameter("KEYFRAME_RATE");

    if(m_global_t>0 && (m_global_t)%kfr==0) {

        vecf m0 = m_motion.back();
        CRigidMotion<float,3> F0inv(m0.Get(3),    // careful with order
                                    m0.Get(4),
                                    m0.Get(5),
                                    m0.Get(0),
                                    m0.Get(1),
                                    m0.Get(2));
        F0inv.Invert();

        vector<pair<vec3f,vec2f> > corrs2i;
        vector<pair<vec2f,vec2f> > corri2i;
        vector<CTracklet*> trackletss2i, trackletsi2i;

        list<shared_ptr<CTracklet> >::iterator it;

        for(it=begin(); it!=end(); it++) {

            vec2f u1 = (*it)->GetLatestLocationAtNativeScale();

            if((*it)->GetStatus() && (*it)->front().HasDescriptor("3DPOINT")) {

                // get access to the scene point
                shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("3DPOINT");
                CDescriptor<vecf>* cdesc = (CDescriptor<vecf>*)desc.get();
                vec3f x(cdesc->Get());

                corrs2i.push_back(pair<vec3f,vec2f>(x,u1));
                trackletss2i.push_back((*it).get());

            }
            else if((*it)->GetStatus() && (*it)->size()>(size_t)kfr) {

                vec2f u0 = (*it)->GetPastLocationAtNativeScale((size_t)kfr);

                corri2i.push_back(pair<vec2f,vec2f>(u0,u1));
                trackletsi2i.push_back((*it).get());

            }

        }


// try 8-point initilialization
//        if(m_global_t/kfr==1) {

//            mat E = CLinearAlgebra::CalibratedNPoint(corri2i,m_cam);
//            mat F = CLinearAlgebra::FactorEssentialMatrix(E);

//            minit(0) = F(0,3);
//            minit(1) = F(1,3);
//            minit(2) = F(2,3);

//            CTransformation<3> T = CTransformation<3>(F);

//            mat R = T.GetLinearPart();

//            vec omega = CRotation<3>::Log(R);

//            minit(3) = omega(0);
//            minit(4) = omega(1);
//            minit(5) = omega(2);
//        }

        cout << "No of correspondences (2d-2d/3d-2d): " << corri2i.size() << " " << corrs2i.size() << endl;

        // init linear solver
        smatf M(0,0);
        CPreconditioner<smatf,vecf,float> precond = CPreconditioner<smatf,vecf,float>(M);
        CIterativeSolver<smatf,vecf,float> solver = CIterativeSolver<smatf,vecf,float>(precond,
                                                                                     m_params->GetIntParameter("CGLS_NITER"),
                                                                                     m_params->GetDoubleParameter("CGLS_EPS"),                                                                                     true);
        // init least-squares problem
        CMagicSfM problem(m_cam,corri2i,corrs2i,F0inv);

        // set up LM method
        CLevenbergMarquardt<smatf,float> lms(problem,solver,m_params->GetDoubleParameter("LM_LAMBDA"));

        // initialize
        vecf& x = problem.Get();

        for(size_t i=0; i<corri2i.size(); i++)
            x(i) = m_params->GetDoubleParameter("INIT_DISTANCE");

        for(size_t i=0; i<6; i++)
            x(corri2i.size()+i) = m0(i);

        // iterate
        vecf r = lms.Iterate(m_params->GetIntParameter("LM_NITER_OUTER"),
                             m_params->GetIntParameter("LM_NITER_INNER"),
                             m_params->GetDoubleParameter("LM_EPS"),
                             false,
                             true);

        // add new points to the map
        double threshold = m_params->GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_DEPTH");

        for(size_t i=0; i<corri2i.size(); i++) {

            if(fabs(r.Get(i))<threshold && fabs(r.Get(corri2i.size()+i))<threshold) {

                // convert depth into world point
                vec3f x0h = m_cam.Normalize(corri2i[i].first);
                vec3f pt = x0h*x(i);
                pt = F0inv.Transform(pt);

                // get pointer to initial feature
                shared_ptr<CAbstractDescriptor> descx(new CDescriptor<vecf>(vecf(pt)));
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
            m1(i) = x(corri2i.size()+i);

        m_motion.push_back(m1);

        // reject outliers
        for(size_t i=0; i<corrs2i.size(); i++) {

            if(fabs(r.Get(2*corri2i.size()+i))>threshold || fabs(r.Get(2*corri2i.size()+corrs2i.size()+i))>threshold)
                trackletss2i[i]->SetStatus(false);

                }


    }

    return 0;

}

bool CMotionTracker::UpdateDescriptors(std::vector<cv::Mat>& pyramid) {

    // get actual pose
    /*vec motion = m_motion.back();
    CRigidMotion<float,3> F(motion.Get(3),           // careful with order
                            motion.Get(4),
                            motion.Get(5),
                            motion.Get(0),
                            motion.Get(1),
                            motion.Get(2));*/

    // smooth image for gradient computation
    vector<Mat> pyrsmooth;
    for(size_t s=0; s<pyramid.size(); s++) {

        Mat imsmooth;
        GaussianBlur(pyramid[s],imsmooth,Size(0,0),m_params->GetDoubleParameter("GRAD_SMOOTH_SIGMA"));

        pyrsmooth.push_back(imsmooth);

    }

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
/*

            // see the feature has an initial depth
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

            }
*/

            // compute brief descriptor
            CBRIEF* briefdesc1 = new CBRIEF(droi);
            shared_ptr<CAbstractDescriptor> brief(briefdesc1);
            brief->Compute(pyrsmooth[s]);
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

#if COMPUTE_ID == 1

            CIdentityDescriptor* tempid = new CIdentityDescriptor(droi,
                                                                  (size_t)m_params->GetIntParameter("NORMALIZE_ID"),
                                                                  (size_t)m_params->GetIntParameter("DESCRIPTOR_HSIZE"));

            tempid->Compute(pyramid[s]);

            shared_ptr<CDescriptor> id(tempid);
            x->AttachDescriptor("ID",id);

#endif

            }

        }

    return 0;

}

bool CMotionTracker::ColorMap(cv::Mat& img) {

    list<shared_ptr<CTracklet> >::iterator it;

    for(it=begin(); it!=end(); it++) {

        if((*it)->GetStatus()) {

            imfeature& x = (*it)->GetLatestState();
            vec2f u0 = x.GetLocationAtNativeScale();

            // add color if not already there
            if((*it)->front().HasDescriptor("3DPOINT") && !(*it)->front().HasDescriptor("COLOR")) {

                size_t row, col;
                row = (size_t)u0(1);
                col = (size_t)u0(0);

                vecf color(3);
                color(0) = (float)img.at<Vec3b>(row,col)[2];
                color(1) = (float)img.at<Vec3b>(row,col)[1];
                color(2) = (float)img.at<Vec3b>(row,col)[0];

                shared_ptr<CAbstractDescriptor> desc(new CDescriptor<vecf>(color));
                (*it)->front().AttachDescriptor("COLOR",desc);

            }

        }

    }

    return 0;

}


CMagicSfM::CMagicSfM(CPinholeCam cam, vector<pair<vec2f,vec2f> >& corri2i, vector<pair<vec3f,vec2f> >& corrs2i, CRigidMotion<float,3> F0inv):
    CLeastSquaresProblem<smatf,float>::CLeastSquaresProblem(2*(corri2i.size()+corrs2i.size()),corri2i.size()+6),
	m_cam(cam),
	m_corri2i(corri2i),
	m_corrs2i(corrs2i),
	m_F0inv(F0inv) {

}

void CMagicSfM::ComputeResidual(vecf& r) {

    // actual transformation from world to second frame
    CRigidMotion<float,3> F1(m_model.Get(m_corri2i.size()+3),           // careful with order
                             m_model.Get(m_corri2i.size()+4),
                             m_model.Get(m_corri2i.size()+5),
                             m_model.Get(m_corri2i.size()),
                             m_model.Get(m_corri2i.size()+1),
                             m_model.Get(m_corri2i.size()+2));

	// projection error for image-to-image correspondences
	for(size_t i=0; i<m_corri2i.size(); i++) {

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

	}

}

void CMagicSfM::ComputeResidualAndJacobian(vecf& r, smatf& J) {

    size_t m = m_corri2i.size();
    size_t n = m_corrs2i.size();

    // actual transformation from world to second frame
    CRigidMotion<float,3> F1(m_model.Get(m_corri2i.size()+3),           // careful with order
                             m_model.Get(m_corri2i.size()+4),
                             m_model.Get(m_corri2i.size()+5),
                             m_model.Get(m_corri2i.size()),
                             m_model.Get(m_corri2i.size()+1),
                             m_model.Get(m_corri2i.size()+2));

    // rotational part of the two vantage points
    matf R0inv = m_F0inv.GetJacobian();
    matf R1 = F1.GetJacobian();

    // derivatives of rotations
    matf DRx(3,3), DRy(3,3), DRz(3,3);
    CDifferentialRotation<float,3>::Rodrigues(m_model.Get(m_corri2i.size()+3),           // careful with order
                                              m_model.Get(m_corri2i.size()+4),
                                              m_model.Get(m_corri2i.size()+5),
                                              nullptr,
                                              DRx.Data().get(),
                                              DRy.Data().get(),
                                              DRz.Data().get());

    // projection error for image-to-image correspondences
    for(size_t i=0; i<m; i++) {

        double wi, wim;
        wi = m_weights.Get(i);
        wim = m_weights.Get(m+i);

        // pixel in img 0 and 1
        vec2f p0, p1;
        p0 = m_corri2i[i].first;
        p1 = m_corri2i[i].second;

        // estimate of the scene point x0
        vec3f x0h = m_cam.Normalize(p0);
        vec3f x0 = x0h*(float)m_model.Get(i);
        x0 = m_F0inv.Transform(x0);       // x0 in world coordinates
        x0h = R0inv*x0h;                  // viewing direction in world coordinates

        // point after rigid transform
        vec3f x1 = F1.Transform(x0);

        // projection of x0 into second image
        vec2f p1p;
        matf Jpi(2,3);
        m_cam.Project(x1,p1p,Jpi);

        // projection error
        vec2f dp = p1p - p1;

        // set residual
        r(i) = wi*dp.Get(0);
        r(m+i) = wim*dp.Get(1);

        // depth derivatives
        vec3f dz = R1*x0h;
        J(i,i) = wi*(Jpi.Get(0,0)*dz.Get(0) + Jpi.Get(0,2)*dz.Get(2));
        J(m+i,i) = wim*(Jpi.Get(1,1)*dz.Get(1) + Jpi.Get(1,2)*dz.Get(2));

        // translational derivative
        J(i,m) = wi*Jpi(0,0); 		J(i,m+1) = wi*Jpi(0,1);	J(i,m+2) = wi*Jpi(0,2);
        J(m+i,m) = wim*Jpi(1,0); 	J(m+i,m+1) = wim*Jpi(1,1); 	J(m+i,m+2) = wim*Jpi(1,2);


        // rotational derivatives.
        vec3f do1 = DRx*x0;
        vec3f do2 = DRy*x0;
        vec3f do3 = DRz*x0;

        J(i,m+3) = wi*(Jpi.Get(0,0)*do1.Get(0) + Jpi.Get(0,2)*do1.Get(2));
        J(m+i,m+3) = wim*(Jpi.Get(1,1)*do1.Get(1) + Jpi.Get(1,2)*do1.Get(2));

        J(i,m+4) = wi*(Jpi.Get(0,0)*do2.Get(0) + Jpi.Get(0,2)*do2.Get(2));
        J(m+i,m+4) = wim*(Jpi.Get(1,1)*do2.Get(1) + Jpi.Get(1,2)*do2.Get(2));

        J(i,m+5) = wi*(Jpi.Get(0,0)*do3.Get(0) + Jpi.Get(0,2)*do3.Get(2));
        J(m+i,m+5) = wim*(Jpi.Get(1,1)*do3.Get(1) + Jpi.Get(1,2)*do3.Get(2));

    }

    // projection error for scene-to-image correspondences
    for(size_t i=0; i<n; i++) {

        double w2mi, w2min;
        w2mi = m_weights.Get(2*m+i);
        w2min = m_weights.Get(2*m+n+i);

        vec3f x0 = m_corrs2i[i].first;		// scene point in world coordinates
        vec2f p1 = m_corrs2i[i].second;		// projection into second image

        // rigid transform from wc to second cam
        vec3f x1 = F1.Transform(x0);

        // projection of x1 into second image
        vec2f p1p;
        matf Jpi(2,3);
        m_cam.Project(x1,p1p,Jpi);

        // projection error
        vec2f dp = p1p - p1;

        // set residual
        r(2*m+i) = w2mi*dp.Get(0);
        r(2*m+n+i) = w2min*dp.Get(1);

        // translational derivative
        J(2*m+i,m) = w2mi*Jpi(0,0);      J(2*m+i,m+1) = w2mi*Jpi(0,1);		J(2*m+i,m+2) = w2mi*Jpi(0,2);
        J(2*m+n+i,m) = w2min*Jpi(1,0); 	J(2*m+n+i,m+1) = w2min*Jpi(1,1);      J(2*m+n+i,m+2) = w2min*Jpi(1,2);

        // rotational derivatives.
        vec3f do1 = DRx*x0;
        vec3f do2 = DRy*x0;
        vec3f do3 = DRz*x0;

        J(2*m+i,m+3) = w2mi*(Jpi.Get(0,0)*do1.Get(0) + Jpi.Get(0,2)*do1.Get(2));
        J(2*m+n+i,m+3) = w2min*(Jpi.Get(1,1)*do1.Get(1) + Jpi.Get(1,2)*do1.Get(2));

        J(2*m+i,m+4) = w2mi*(Jpi.Get(0,0)*do2.Get(0) + Jpi.Get(0,2)*do2.Get(2));
        J(2*m+n+i,m+4) = w2min*(Jpi.Get(1,1)*do2.Get(1) + Jpi.Get(1,2)*do2.Get(2));

        J(2*m+i,m+5) = w2mi*(Jpi.Get(0,0)*do3.Get(0) + Jpi.Get(0,2)*do3.Get(2));
        J(2*m+n+i,m+5) = w2min*(Jpi.Get(1,1)*do3.Get(1) + Jpi.Get(1,2)*do3.Get(2));

    }

}

vecf CMagicSfM::ComputeDispersion(const vecf& r) {

    vecf sigma(r.NElems());

    size_t m = m_corri2i.size();
    size_t n = m_corrs2i.size();

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
