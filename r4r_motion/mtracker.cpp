/*
 * mtracker.cpp
 *
 *  Created on: Jul 26, 2012
 *      Author: jbalzer
 */




#include "mtracker.h"
#include "lm.h"
#include "basic.h"
#include "utils.h"
#include "trafo.h"
#include "interp.h"


#include <fstream>

using namespace cv;

namespace R4R {


CMotionTracker::CMotionTracker(CParameters params, CCam cam):
		CSimpleTracker(params),
		m_cam(cam),
		m_motion() {

		m_motion.push_back(vec(6));

}

map<CTracklet*,vec> CMotionTracker::GetMap() {

	map<CTracklet*,vec> pts;
	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			vec ptcolor(6);

			// check whether we have 3d data for feature in first frame
            if((*it)->front().HasDescriptor("3DPOINT")) {

				// get 3d point
                shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("3DPOINT");
                CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
				vec pt = cdesc->Get();

				// copy data
				ptcolor(0) = pt(0);
				ptcolor(1) = pt(1);
				ptcolor(2) = pt(2);

			}

			// check whether the point is colored
            if((*it)->front().HasDescriptor("COLOR")) {

				// color
                shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("COLOR");
                CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
				vec color = cdesc->Get();

				// copy data
				ptcolor(3) = color(0);
				ptcolor(4) = color(1);
				ptcolor(5) = color(2);

			}

			// only attach if it has at least a point
            if((*it)->front().HasDescriptor("3DPOINT"))
				pts.insert(pair<CTracklet*,vec>((*it).get(),ptcolor));

		}

	}

	return pts;

}

bool CMotionTracker::SaveMotion(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file " << filename << "." << endl;
		return 1;

	 }

	list<vec>::iterator it;

	for(it=m_motion.begin(); it!=m_motion.end(); it++) {

		vec m = (*it);
		m.Transpose();
		out << m << endl;

	}

	out.close();

	return 0;

}

bool CMotionTracker::SaveMap(const char* filename) {

	map<CTracklet*,vec> pts = GetMap();
	map<CTracklet*,vec>::iterator it;

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file " << filename << "." << endl;
		return 1;

	 }

	out << "ply" << endl;
	out <<  "format ascii 1.0" << endl;
	out <<  "comment" << endl;
	out << "element vertex " << pts.size() << endl;
	out << "property float32 x" << endl;
	out << "property float32 y" << endl;
	out << "property float32 z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "end_header" << endl;

	for(it=pts.begin(); it!=pts.end(); it++) {

		vec cpt = it->second;

		out << cpt.Get(0) << " " << cpt.Get(1) << " " << cpt.Get(2) << " " << (unsigned int)cpt.Get(3) << " " <<  (unsigned int)cpt.Get(4) << " " <<  (unsigned int)cpt.Get(5) << endl;

	}

	out.close();

	return 0;

}

vec CMotionTracker::Triangulate(vec o0, vec d0, vec o1, vec d1, double eps) {

	vec result;

	vec z = o0 - o1;

	double ab, aa, bb, az, bz;
	ab = vec::InnerProduct(d0,d1);
	aa = vec::InnerProduct(d0,d0);
	bb = vec::InnerProduct(d1,d1);
	az = vec::InnerProduct(d0,z);
	bz = vec::InnerProduct(d1,z);

	double lambda = (ab*bz-bb*az)/(aa*bb-ab*ab);
	double mu = -(ab*az-aa*bz)/(aa*bb-ab*ab);

	vec x0s = o0 + d0*lambda;
	vec x1s = o1 + d1*mu;
	double dist = (x1s - x0s).Norm2();

	if(lambda>0 && mu>0 && dist<eps)
		result = (x0s + x1s)*0.5;

	return result;

}


vec CMotionTracker::EpipolarSearch(vec& ub, vec& d, cv::Mat& img0, cv::Mat& img1, size_t hsize, size_t nsteps, double eps) {

	size_t w = 2*hsize + 1;
	size_t n = w*w;

	vec r(n);
	vec J(n);

	double rold, rnew, lambda;

	rnew = 0;
	lambda = 0;

	vec u1;

	for(size_t k=0; k<nsteps; k++) {

		// predict new target point location
		u1 = ub + d*lambda;

		for(int i=-(int)hsize; i<=(int)hsize; i++) {

			for(int j=-(int)hsize; j<=(int)hsize; j++) {

				// compute row index
				int row = w*(i+(int)hsize)+(int)hsize+j;

				// interpolated intensities and gradients
				double I0, I1, I1u, I1v;
				I0 = CImageInterpolation::Bilinear(img0,ub(0)+i,ub(1)+j);			// FIXME: this does not have to be recomputed everytime!!!!
				I1 = CImageInterpolation::Bilinear(img1,u1(0)+i,u1(1)+j);
				I1u = CImageInterpolation::Gradient(img1,u1(0)+i,u1(1)+j,true);
				I1v = CImageInterpolation::Gradient(img1,u1(0)+i,u1(1)+j,false);

				// residual
				r(row) = I1 - I0;
				J(row) = I1u*d(0) + I1v*d(1);

			}

		}

		// save old residual norm and compute new one
		rold = rnew;
		rnew = r.Norm2();

		// break if there is no more change
		if(fabs(rnew-rold)<eps)
			break;

		// GN step
		double dlambda = -vec::InnerProduct(J,r)/vec::InnerProduct(J,J);
		lambda = lambda + dlambda;

	}

	return u1;

}


bool CMotionTracker::Update(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

	// do LK tracking at every time instance, FIXME: how to do descriptor aggregation with 3d info
	CSimpleTracker::Update(pyramid0,pyramid1);

	size_t kfr = (size_t)m_params.GetIntParameter("KEYFRAME_RATE");

	if(m_global_t>0 && (m_global_t)%kfr==0) {

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//		stringstream no;
//		no.fill('0');
//		no.width(4);
//		no << m_global_t;
//
//		string prefix0("/home/jbalzer/Dump/corri2i");
//		string prefix1("/home/jbalzer/Dump/corrs2i");
//		stringstream filename0, filename1;
//		filename0 << prefix0 << no.str() << ".txt";
//		filename1 << prefix1 << no.str() << ".txt";
//
//		mat corr0, corr1;

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		vec m0 = m_motion.back();
		mat F0 = CLinearAlgebra::CreateTransformationMatrix(m0);
		mat F0inv = CLinearAlgebra::InvertTransformation(F0);

		vector<pair<vec,vec> > corrs2i, corri2i;
		vector<CTracklet*> trackletss2i, trackletsi2i;

		list<shared_ptr<CTracklet> >::iterator it;

		for(size_t s=0; s<size(); s++) {

			for(it=at(s).begin(); it!=at(s).end(); it++) {

				vec u1 = (*it)->GetLatestLocationAtNativeScale();

                if((*it)->GetStatus() && (*it)->front().HasDescriptor("3DPOINT")) {

					// get access to the scene point
                    shared_ptr<CAbstractDescriptor> desc = (*it)->front().GetDescriptor("3DPOINT");
                    CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
					vec x = cdesc->Get();

					corrs2i.push_back(pair<vec,vec>(x,u1));
					trackletss2i.push_back((*it).get());

				}
				else if((*it)->GetStatus() && (*it)->size()>(size_t)kfr) {

					vec u0 = (*it)->GetPastLocationAtNativeScale((size_t)kfr);

					corri2i.push_back(pair<vec,vec>(u0,u1));
					trackletsi2i.push_back((*it).get());

				}

			}

		}

		/*if(m_global_t/kfr==1) {

			mat E = CLinearAlgebra::CalibratedNPoint(corri2i,m_cam);
			mat F = CLinearAlgebra::FactorEssentialMatrix(E);

			minit(0) = F(0,3);
			minit(1) = F(1,3);
			minit(2) = F(2,3);

			CTransformation<3> T = CTransformation<3>(F);

			mat R = T.GetLinearPart();

			vec omega = CRotation<3>::Log(R);

			minit(3) = omega(0);
			minit(4) = omega(1);
			minit(5) = omega(2);
		}*/

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////


//		corr0 = mat(corri2i.size(),9); // correspondence, two residuals, new 3d point
//		corr1 = mat(corrs2i.size(),6); // correspondence, two residuals
//
//		for(size_t i=0; i<corri2i.size(); i++) {
//
//			corr0(i,0) = corri2i[i].first.Get(0);
//			corr0(i,1) = corri2i[i].first.Get(1);
//			corr0(i,2) = corri2i[i].second.Get(0);
//			corr0(i,3) = corri2i[i].second.Get(1);
//
//
//		}
//
//		for(size_t i=0; i<corrs2i.size(); i++) {
//
//			corr1(i,0) = corrs2i[i].first.Get(0);
//			corr1(i,1) = corrs2i[i].first.Get(1);
//			corr1(i,2) = corrs2i[i].second.Get(0);
//			corr1(i,3) = corrs2i[i].second.Get(1);
//
//		}


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////

		cout << "No of correspondences (2d-2d/3d-2d): " << corri2i.size() << " " << corrs2i.size() << endl;

		// init linear solver
		smat M(0,0);
		CPreconditioner<smat,vec,double> precond = CPreconditioner<smat,vec,double>(M);
		CIterativeSolver<smat,vec,double> solver = CIterativeSolver<smat,vec,double>(precond,
																					 m_params.GetIntParameter("CGLS_NITER"),
																					 m_params.GetDoubleParameter("CGLS_EPS"),
																					 true);

		// init least-squares problem
		CMagicSfM problem(m_cam,corri2i,corrs2i,F0inv);

		// set up LM method, FIXME: get params from file
		CLevenbergMarquardt<smat> lms(problem,solver,m_params.GetDoubleParameter("LM_LAMBDA"));

		// initialize
		vec& x = problem.Get();

		for(size_t i=0; i<corri2i.size(); i++)
			x(i) = m_params.GetDoubleParameter("INIT_DISTANCE");

		for(size_t i=0; i<6; i++)
			x(corri2i.size()+i) = m0(i);

		// iterate
		vec r = lms.Iterate(m_params.GetIntParameter("LM_NITER_OUTER"),
							m_params.GetIntParameter("LM_NITER_INNER"),
							m_params.GetDoubleParameter("LM_EPS"),
							false,
							true);

		// add new points to the map
		double threshold = m_params.GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_DEPTH");

		for(size_t i=0; i<corri2i.size(); i++) {

			if(fabs(r.Get(i))<threshold && fabs(r.Get(corri2i.size()+i))<threshold) {

				// convert depth into world point
				vec x0h = m_cam.NormalizeLocal(corri2i[i].first);
				vec pt = x0h*x(i);
				pt = CLinearAlgebra::TransformPoint(F0inv,pt);		// x0 in world coordinates

				///////////////////////////////////////////////////////////////////////////////////////////////////////

				//corr0(i,4) = r.Get(i);
				//corr0(i,5) = r.Get(corri2i.size()+i);
				//corr0(i,6) = pt(0);
				//corr0(i,7) = pt(1);
				//corr0(i,8) = pt(2);

				///////////////////////////////////////////////////////////////////////////////////////////////////////

				// get pointer to initial feature
                CFeature feat = trackletsi2i[i]->front();
                shared_ptr<CAbstractDescriptor> descx(new CDescriptor<vec>(pt));
                feat.AttachDescriptor("3DPOINT",descx);

				// also store initial depths for descriptor canonization
				/*vec z(1);
				z(0) = x(i);
                shared_ptr<CDescriptor> descz(new CDescriptor<vec>(z));
				feat->AttachDescriptor("INITDEPTH",descz);*/

			}
			else
				trackletsi2i[i]->SetStatus(false);

		}

		// process motion
		threshold = m_params.GetDoubleParameter("OUTLIER_REJECTION_THRESHOLD_MOTION");

		vec m1(6);
		for(size_t i=0; i<6; i++)
			m1(i) = x(corri2i.size()+i);


		m_motion.push_back(m1);


		// reject outliers
		for(size_t i=0; i<corrs2i.size(); i++) {

			if(fabs(r.Get(2*corri2i.size()+i))>threshold || fabs(r.Get(2*corri2i.size()+corrs2i.size()+i))>threshold)
				trackletss2i[i]->SetStatus(false);

			///////////////////////////////////////////////////////////////////////////////////////////////////////

			//corr1(i,4) = r.Get(corri2i.size()+i);
			//corr1(i,5) = r.Get(2*corri2i.size()+corrs2i.size()+i);

			///////////////////////////////////////////////////////////////////////////////////////////////////////

		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//corr0.SaveToFile(filename0.str().c_str());
		//corr1.SaveToFile(filename1.str().c_str());

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	}


	return 0;

}

bool CMotionTracker::UpdateDescriptors(std::vector<cv::Mat>& pyramid) {

	list<shared_ptr<CTracklet> >::iterator it;

	vec motion = m_motion.back();
	mat F = CLinearAlgebra::CreateTransformationMatrix(motion);

	for(size_t s=0; s<size(); s++) {

		Mat imsmooth;
		GaussianBlur(pyramid[s],imsmooth,Size(0,0),m_params.GetDoubleParameter("GRAD_SMOOTH_SIGMA"));

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

                CFeature& x = (*it)->GetLatestState();
                vec u0 = x.GetLocation();

				// create new feature
				CRectangle<double> droi(u0.Get(0),
										u0.Get(1),
									    m_params.GetIntParameter("DESCRIPTOR_HSIZE"),
									    m_params.GetIntParameter("DESCRIPTOR_HSIZE"));

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

#if COMPUTE_ID == 1
				CIdentityDescriptor* tempid = new CIdentityDescriptor(droi,
																	  (size_t)m_params.GetIntParameter("NORMALIZE_ID"),
																	  (size_t)m_params.GetIntParameter("DESCRIPTOR_HSIZE"));

				tempid->Compute(pyramid[s]);

				shared_ptr<CDescriptor> id(tempid);
				x->AttachDescriptor("ID",id);
#endif

				}

			}

		}

	return 0;

}

bool CMotionTracker::ColorMap(cv::Mat& img) {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

                CFeature& x = (*it)->GetLatestState();
                vec u0 = x.GetLocationAtNativeScale();

				// add color if not already there
                if((*it)->front().HasDescriptor("3DPOINT") && !(*it)->front().HasDescriptor("COLOR")) {

					size_t row, col;
					row = (size_t)u0(1);
					col = (size_t)u0(0);

					vec color(3);
					color(0) = (double)img.at<Vec3b>(row,col)[2];
					color(1) = (double)img.at<Vec3b>(row,col)[1];
					color(2) = (double)img.at<Vec3b>(row,col)[0];

                    CFeature feat = (*it)->front();
                    shared_ptr<CAbstractDescriptor> desc(new CDescriptor<vec>(color));
                    feat.AttachDescriptor("COLOR",desc);

				}

			}

		}

	}

	return 0;

}


CMagicSfM::CMagicSfM(CCam cam, vector<pair<vec,vec> >& corri2i, vector<pair<vec,vec> >& corrs2i, mat F0inv):
	CLeastSquaresProblem<smat>::CLeastSquaresProblem(2*(corri2i.size()+corrs2i.size()),corri2i.size()+6),
	m_cam(cam),
	m_corri2i(corri2i),
	m_corrs2i(corrs2i),
	m_F0inv(F0inv) {

}

void CMagicSfM::ComputeResidual(vec& r) {

	// get actual transformation, six motion parameters are stored at the end
	vec t(3), omega(3);
	for(size_t i=0; i<3; i++) {

		t(i) = m_model(m_corri2i.size()+i);
		omega(i) = m_model(m_corri2i.size()+3+i);

	}

	mat R = CLinearAlgebra::Rodrigues(omega);

	// projection error for image-to-image correspondences
#pragma omp parallel for
	for(size_t i=0; i<m_corri2i.size(); i++) {

		// pixel in img 0 and 1
		vec u0, u1;
		u0 = m_corri2i[i].first;
		u1 = m_corri2i[i].second;

		// estimate of the scene point x0 in camera system of first frame
		vec x0h = m_cam.NormalizeLocal(u0);
		vec x0 = x0h*m_model(i);
		x0 = CLinearAlgebra::TransformPoint(m_F0inv,x0);

		// point transformed from world to local camera system of second frame
		vec x1 = R*x0 + t;

		// projection of x0 into second image
		vec u1p(2);
		mat Jpi(2,3);
		m_cam.ProjectLocal(x1,u1p,Jpi);

		// projection error
		vec du = u1p - u1;

		// set residual
		r(i) = du(0);
		r(m_corri2i.size()+i) = du(1);

	}

	// projection error for scene-to-image correspondences
#pragma omp parallel for
	for(size_t i=0; i<m_corrs2i.size(); i++) {

		vec x0, u1;
		x0 = m_corrs2i[i].first;		// scene point in world coordinates
		u1 = m_corrs2i[i].second;		// projection into second image

		// rigid transform from wc to second cam
		vec x1 = R*x0 + t;

		// projection of x1 into second image
		vec u1p(2);
		mat Jpi(2,3);
		m_cam.ProjectLocal(x1,u1p,Jpi);

		// projection error
		vec du = u1p - u1;

		// set residual
		r(2*m_corri2i.size()+i) = du(0);
		r(2*m_corri2i.size()+m_corrs2i.size()+i) = du(1);

	}

}

void CMagicSfM::ComputeResidualAndJacobian(vec& r, smat& J) {

	size_t m = m_corri2i.size();
	size_t n = m_corrs2i.size();

	// get actual transformation, six motion parameters are stored at the end
	vec t(3), omega(3);
	for(size_t i=0; i<3; i++) {

		t(i) = m_model(m+i);
		omega(i) = m_model(m+3+i);

	}

	mat R(3,3), DRx(3,3), DRy(3,3), DRz(3,3);
	CLinearAlgebra::Rodrigues(omega,R,DRx,DRy,DRz);

	// projection error for image-to-image correspondences
//#pragma omp parallel for
	for(size_t i=0; i<m_corri2i.size(); i++) {

		// pixel in img 0 and 1
		vec u0, u1;
		u0 = m_corri2i[i].first;
		u1 = m_corri2i[i].second;

		// estimate of the scene point x0, FIXME: this is in local coordinates, just work at initialization
		vec x0h = m_cam.NormalizeLocal(u0);
		vec x0 = x0h*m_model(i);
		x0 = CLinearAlgebra::TransformPoint(m_F0inv,x0);		// x0 in world coordinates
		x0h = CLinearAlgebra::TransformDirection(m_F0inv,x0h);	// viewing direction in world coordinates

		// point after rigid transform
		vec x1 = R*x0 + t;

		// projection of x0 into second image
		vec u1p(2);
		mat Jpi(2,3);
		m_cam.ProjectLocal(x1,u1p,Jpi);

		// projection error
		vec du = u1p - u1;

		// set residual
		r(i) = du(0);
		r(m+i) = du(1);

		// depth derivatives
		vec dz = Jpi*(R*x0h);
		J(i,i) = dz(0);
		J(m+i,i) = dz(1);

		// translational derivative
		J(i,m) = Jpi(0,0); 		J(i,m+1) = Jpi(0,1);	J(i,m+2) = Jpi(0,2);
		J(m+i,m) = Jpi(1,0); 	J(m+i,m+1) = Jpi(1,1); 	J(m+i,m+2) = Jpi(1,2);

		// rotational derivatives.
		vec do1 = Jpi*(DRx*x0);
		vec do2 = Jpi*(DRy*x0);
		vec do3 = Jpi*(DRz*x0);

		J(i,m+3) = do1(0); 		J(i,m+4) = do2(0); 		J(i,m+5) = do3(0);
		J(m+i,m+3) = do1(1); 	J(m+i,m+4) = do2(1);	J(m+i,m+5) = do3(1);

	}


	// projection error for scene-to-image correspondences
//#pragma omp parallel for
	for(size_t i=0; i<n; i++) {

		vec x0, u1;
		x0 = m_corrs2i[i].first;		// scene point in world coordinates
		u1 = m_corrs2i[i].second;		// projection into second image

		// rigid transform from wc to second cam
		vec x1 = R*x0 + t;

		// projection of x1 into second image
		vec u1p(2);
		mat Jpi(2,3);
		m_cam.ProjectLocal(x1,u1p,Jpi);

		// projection error
		vec du = u1p - u1;

		// set residual
		r(2*m+i) = du(0);
		r(2*m+n+i) = du(1);

		// translational derivative
		J(2*m+i,m) = Jpi(0,0); 	J(2*m+i,m+1) = Jpi(0,1);		J(2*m+i,m+2) = Jpi(0,2);
		J(2*m+n+i,m) = Jpi(1,0); 	J(2*m+n+i,m+1) = Jpi(1,1); 	J(2*m+n+i,m+2) = Jpi(1,2);

		// rotational derivatives.
		vec do1 = Jpi*(DRx*x0);
		vec do2 = Jpi*(DRy*x0);
		vec do3 = Jpi*(DRz*x0);

		J(2*m+i,m+3) = do1(0); 		J(2*m+i,m+4) = do2(0);		J(2*m+i,m+5) = do3(0);
		J(2*m+n+i,m+3) = do1(1);	J(2*m+n+i,m+4) = do2(1); 	J(2*m+n+i,m+5) = do3(1);

	}

}


vec CMagicSfM::ComputeDispersion(const vec& r) {

	vec sigma(r.NElems());

	size_t m = m_corri2i.size();
	size_t n = m_corrs2i.size();

	// get pointer to the data
	double* rdata = r.Data();

	// copy into subarrays (one for each dimension and part of the functional)
	vec ri2iu(m,&rdata[0]);
	vec ri2iv(m,&rdata[m]);
	vec rs2iu(n,&rdata[2*m]);
	vec rs2iv(n,&rdata[2*m+n]);

	// get a robust estimate of the mean
/*
	vec medi2i(2), meds2i(2);
	medi2i(0) = ri2iu.Median();
	medi2i(1) = ri2iv.Median();
	meds2i(0) = rs2iu.Median();
	meds2i(1) = rs2iv.Median();
*/

	double si2iu = 1.0/(1.4826*CLinearAlgebra::MedianAbsoluteDeviation(ri2iu));
	double si2iv = 1.0/(1.4826*CLinearAlgebra::MedianAbsoluteDeviation(ri2iv));
	double ss2iu = 1.0/(1.4826*CLinearAlgebra::MedianAbsoluteDeviation(rs2iu));
	double ss2iv = 1.0/(1.4826*CLinearAlgebra::MedianAbsoluteDeviation(rs2iv));

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





CSfMTrackerUpdate::CSfMTrackerUpdate(CCam cam, vector<CTracklet*> map, cv::Mat& img0, cv::Mat& img1, size_t hsize):
	CLeastSquaresProblem(map.size()*(2*hsize+1)*(2*hsize+1),6),
	m_cam(cam),
	m_map(map),
	m_img0(img0),
	m_img1(img1),
	m_hsize(hsize) {


}

void CSfMTrackerUpdate::ComputeResidualAndJacobian(vec& r, mat& J) {

	// some sizes
	size_t n = m_map.size();					// no of anchor pixels
	size_t w = 2*m_hsize+1; 					// width of neighborhood
	size_t m = w*w;								// no of pixels around each anchor

	// actual transformation
	vec t(3), omega(3);
	for(size_t i=0; i<3; i++) {

		t(i) = m_model(i);
		omega(i) = m_model(3+i);

	}

	mat R(3,3), DRx(3,3), DRy(3,3), DRz(3,3);
	CLinearAlgebra::Rodrigues(omega,R,DRx,DRy,DRz);

	// projection errors and their derivatives
	for(size_t i=0; i<n; i++) {

		// anchor point in first frame
		vec u0 = m_map[i]->GetLatestLocation();

		// world point
        shared_ptr<CAbstractDescriptor> desc = m_map[i]->front().GetDescriptor("3DPOINT");
        CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
		vec x0 = cdesc->Get();

		// rigid transform from wc to second cam
		vec x1 = R*x0 + t;

		// projection of x1 into second image
		vec u1 = m_cam.ProjectLocal(x1);

		// init Jacobian w.r.t to local coordinates of projection into second frame
		mat Jpi(2,3);
		Jpi(0,0) = m_cam.m_f[0]/x1(2);
		Jpi(1,1) = m_cam.m_f[1]/x1(2);

		// every pixel in neighborhood of anchor contributes to energy
		for(int j=-(int)m_hsize; j<=(int)m_hsize; j++) {

			for(int k=-(int)m_hsize; k<=(int)m_hsize; k++) {

				// compute row index
				int row = m*i+w*(j+(int)m_hsize)+(int)m_hsize+k;

				// interpolated intensities and gradients
				double I0, I1, I1u, I1v;
				I0 = CImageInterpolation::Bilinear(m_img0,u0(0)+j,u0(1)+k);
				I1 = CImageInterpolation::Bilinear(m_img1,u1(0)+j,u1(1)+k);
				I1u = CImageInterpolation::Gradient(m_img1,u1(0)+j,u1(1)+k,true);
				I1v = CImageInterpolation::Gradient(m_img1,u1(0)+j,u1(1)+k,false);

				// residual
				r(row) = I1 - I0;

				// update projection Jacobian
				Jpi(0,2) = -(u1(0)+j-m_cam.m_c[0])/x1(2);
				Jpi(1,2) = -(u1(1)+k-m_cam.m_c[1])/x1(2);

				// translational part
				J(row,0) = I1u*Jpi(0,0) + I1v*Jpi(1,0);
				J(row,1) = I1u*Jpi(0,1) + I1v*Jpi(1,1);
				J(row,2) = I1u*Jpi(0,2) + I1v*Jpi(1,2);

				// rotational part
				vec do1 = Jpi*(DRx*x0);
				vec do2 = Jpi*(DRy*x0);
				vec do3 = Jpi*(DRz*x0);
				J(row,3) = I1u*do1(0) + I1v*do1(1);
				J(row,4) = I1u*do2(0) + I1v*do2(1);
				J(row,5) = I1u*do3(0) + I1v*do3(1);

			}

		}

	}

}

void CSfMTrackerUpdate::ComputeResidual(vec& r) {
	// some sizes
	size_t n = m_map.size();					// no of anchor pixels
	size_t w = 2*m_hsize+1; 					// width of neighborhood
	size_t m = w*w;								// no of pixels around each anchor

	// actual transformation
	vec t(3), omega(3);
	for(size_t i=0; i<3; i++) {

		t(i) = m_model(i);
		omega(i) = m_model(3+i);

	}

	mat R(3,3), DRx(3,3), DRy(3,3), DRz(3,3);
	CLinearAlgebra::Rodrigues(omega,R,DRx,DRy,DRz);

	// projection errors and their derivatives
#pragma omp parallel for
	for(size_t i=0; i<n; i++) {

		// anchor point in first frame
		vec u0 = m_map[i]->GetLatestLocation();

		// world point
        shared_ptr<CAbstractDescriptor> desc = m_map[i]->front().GetDescriptor("3DPOINT");
        CDescriptor<vec>* cdesc = (CDescriptor<vec>*)desc.get();
		vec x0 = cdesc->Get();

		// rigid transform from wc to second cam
		vec x1 = R*x0 + t;

		// projection of x1 into second image
		vec u1 = m_cam.ProjectLocal(x1);

		// every pixel in neighborhood of anchor contributes to energy
		for(int j=-(int)m_hsize; j<=(int)m_hsize; j++) {

			for(int k=-(int)m_hsize; k<=(int)m_hsize; k++) {

				// compute row index
				int row = m*i+w*(j+(int)m_hsize)+(int)m_hsize+k;

				// interpolated intensities and gradients
				double I0, I1;
				I0 = CImageInterpolation::Bilinear(m_img0,u0(0)+j,u0(1)+k);	// FIXME: do this outside
				I1 = CImageInterpolation::Bilinear(m_img1,u1(0)+j,u1(1)+k);

				// residual
				r(row) = I1 - I0;


			}

		}

	}

}


}
