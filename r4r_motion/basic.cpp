/*
 * iddesc.cpp
 *
 *  Created on: May 9, 2012
 *      Author: jbalzer
 */




#include "basic.h"
#include "rect.h"
#include "interp.h"

using namespace std;
using namespace cv;

namespace R4R {

CIdentityDescriptor::CIdentityDescriptor(CRectangle<double> roi, size_t method, size_t hsize):
    CNeighborhoodDescriptor(matf(2*hsize+1,2*hsize+1),roi),
    m_method(method) {}

bool CIdentityDescriptor::Compute(cv::Mat& img) {

	double dx = 2.0/(double)(m_container.NCols()-1);
	double dy = 2.0/(double)(m_container.NRows()-1);

	for(size_t i=0; i<m_container.NRows(); i++) {

		for(size_t j=0; j<m_container.NCols(); j++) {

			vec pt = m_roi.TransformFrom(-1.0+j*dx,-1.0+i*dy);

			m_container(i,j) = CImageInterpolation::Bilinear(img,pt(0),pt(1));

		}

	}

	switch(m_method) {

	case 1:
		NormalizeWeber();
		break;
	case 2:
		Normalize();
		break;
	case 3:
		Center();
		break;
	default:
		break;

	}

	return 0;

}

void CIdentityDescriptor::Center() {

	double sum = m_container.Sum();

	if(sum==0)
		return;

	double mean = sum/m_container.NElems();

	m_container = m_container - mean;


}

void CIdentityDescriptor::NormalizeWeber() {

	double sum = m_container.Sum();

	double mean = sum/m_container.NElems();

	if(mean==0)
		return;

	m_container.Scale(1.0/mean);

}


void CIdentityDescriptor::Normalize() {

	Center();

	m_container.Normalize();

}

double CIdentityDescriptor::Distance(CDescriptor<matf>& desc) const {

    return matf::InnerProduct(this->m_container,desc.Get());

}

CIdentityGradientDescriptor::CIdentityGradientDescriptor(CRectangle<double> roi, double alpha, size_t method, size_t hsize):
    CNeighborhoodDescriptor(mat(4*hsize+2,2*hsize+1),roi),
	m_alpha(alpha),
	m_method(method) {

}

bool CIdentityGradientDescriptor::Compute(cv::Mat& img) {

	size_t h = m_container.NRows()/2;
	double dx = 2.0/(double)(m_container.NCols()-1);
	double dy = 2.0/(double)(h-1);

	for(size_t i=0; i<h; i++) {

		for(size_t j=0; j<m_container.NCols(); j++) {

			vec pt = m_roi.TransformFrom(-1.0+j*dx,-1.0+i*dy);

			m_container(i,j) = CImageInterpolation::Gradient(img,pt(0),pt(1),0);
			m_container(h+i,j) = CImageInterpolation::Gradient(img,pt(0),pt(1),1);

		}

	}

	if(m_method>0)
		Normalize();

	return 0;

}

void CIdentityGradientDescriptor::Normalize() {

	size_t h = m_container.NRows()/2;

	for(size_t i=0; i<h; i++) {

		for(size_t j=0; j<m_container.NCols(); j++) {

			double gnorm = sqrt(m_container.Get(i,j)*m_container.Get(i,j) + m_container.Get(h+i,j)*m_container.Get(h+i,j));
			double weight = WeightingFunction(gnorm);

			m_container(i,j) *= weight;
			m_container(h+i,j) *= weight;

		}

	}

}

double CIdentityGradientDescriptor::Distance(CDescriptor<mat>& desc) const {

	size_t h = m_container.NRows()/2;

	mat& container = desc.Get();
	mat sp(h,m_container.NCols());

	for(size_t i=0; i<h; i++) {

		for(size_t j=0; j<m_container.NCols(); j++) {

			sp(i,j) = m_container.Get(i,j)*container.Get(h+i,j) - m_container.Get(h+i,j)*container.Get(i,j);

		}

	}

	sp = sp - sp.Mean();

	return sp.Norm1();

}

double CIdentityGradientDescriptor::WeightingFunction(double gnorm) {

	double w;

	switch(m_method) {

	case 1:

		if(gnorm<m_alpha)
			w = 0;
		else
			w = 1.0/gnorm;

		break;

	case 2:

		w = (2*m_alpha*gnorm)/(m_alpha*gnorm*gnorm + 1.0);

		break;

	default:

		w = 1;

		break;

	}

	return w;

}


CCurvatureDescriptor::CCurvatureDescriptor(CRectangle<double> roi, double alpha, size_t method, size_t hsize):
    CIdentityGradientDescriptor(roi,alpha,method,hsize),
    m_kappa(2*hsize+1,2*hsize+1) {}

bool CCurvatureDescriptor::Compute(cv::Mat& img) {

	// compute gradient field
	CIdentityGradientDescriptor::Compute(img);

	// normalize
	Normalize();

	size_t h = m_container.NRows()/2;

	// compute divergence by central differences
	for(size_t i=1; i<h-1; i++) {

		for(size_t j=1; j<m_container.NCols()-1; j++) {

			double dvudu = 0.5*(m_container(i,j+1) - m_container(i,j-1));
			double dvvdv = 0.5*(m_container(h+i+1,j) - m_container(h+i-1,j));
            m_kappa(i-1,j-1) = 0.5*(dvudu + dvvdv);

		}

	}

	return 0;

}

double CCurvatureDescriptor::Distance(const CCurvatureDescriptor& desc) const {

    mat diff = m_kappa - desc.m_kappa;

    return diff.Norm2();

}



double CBRIEF::m_x[2*BITSET_LENGTH] = {};

double CBRIEF::m_y[2*BITSET_LENGTH] = {};


CBRIEF::CBRIEF(CRectangle<double> roi):
    CNeighborhoodDescriptor(CDenseVector<bool>(BITSET_LENGTH),roi)
    {}


void CBRIEF::GenerateSamplePoints() {

    // draw samples on the set [-1,1]x[-1,1]
    mat sigma(2,2);
    sigma.Eye();
    sigma.Scale(2.0/5.0);

    // loop through number of samples
    for(size_t i=0; i<BITSET_LENGTH; i++) {

        // loop through number of points per test
        for(size_t j=0; j<2; j++) {

            // Box-Muller transform
            double x, y;
            x = rand()/double(RAND_MAX);
            y = rand()/double(RAND_MAX);

            double r = sqrt(-2*log(x));

            vec z(2);
            z(0) = r*cos(2*M_PI*y);
            z(1) = r*sin(2*M_PI*y);

            vec p = sigma*z;

            switch(j) {

            case 0:

                CBRIEF::m_x[i] = p(0);
                CBRIEF::m_x[BITSET_LENGTH+i] = p(1);

                break;

            case 1:

                CBRIEF::m_y[i] = p(0);
                CBRIEF::m_y[BITSET_LENGTH+i] = p(1);

                break;

            }

        }

    }

}

bool CBRIEF::Compute(cv::Mat& img) {

    // smoothing
    //Mat img_smooth;
    //cv::GaussianBlur(img,img_smooth,cv::Size(7,7),2);

    for(size_t i=0; i<m_container.NElems(); i++) {

        // get sample points from static member
        vec x(2), y(2);
        x(0) = CBRIEF::m_x[i];
        x(1) = CBRIEF::m_x[BITSET_LENGTH+i];
        y(0) = CBRIEF::m_y[i];
        y(1) = CBRIEF::m_y[BITSET_LENGTH+i];

        // transform to image coordinates
        x = m_roi.TransformFrom(x);
        y = m_roi.TransformFrom(y);

        // interpolate image
        double Ix, Iy;
        Ix = CImageInterpolation::Bilinear(img,x(0),x(1));
        Iy = CImageInterpolation::Bilinear(img,y(0),y(1));

        // perform binary test
        m_container(i) = (Iy>Ix);

    }

    return 0;

}

double CBRIEF::Distance(const CBRIEF& desc) const {

    CDenseVector<bool> diff = m_container - desc.m_container;

    return diff.HammingNorm();

}

/*
CFBDDescriptor::CFBDDescriptor(CRectangle<double> roi, size_t length):
	CNeighborhoodDescriptor(roi),
	m_length(length) {

}

bool CFBDDescriptor::Compute(cv::Mat& img) {

	// FIXME: generate sample points in static function to reduce computation

	//double r =

	m_container = vec(m_length);

	double r = min(m_roi.Width(),m_roi.Height());



	return 0;

}*/


}
