/*
 * cam.cpp
 *
 *  Created on: May 30, 2012
 *      Author: jbalzer
 */


#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#include "cam.h"
#include "rutils.h"
#include "factor.h"

using namespace std;

namespace R4R {

CCam::CCam() {

	m_size[0] = 0;
	m_size[1] = 0;
	m_f[0] = 0;
	m_f[1] = 0;
	m_c[0] = 0;
	m_c[1] = 0;

	for(size_t i=0;i<5;i++)
		m_k[i]=0;

	m_alpha = 0;

	m_F = mat(4,4);
	m_F.Eye();

	m_Finv = m_F;

}

CCam::CCam(size_t w, size_t h) {

	m_size[0] = w;
	m_size[1] = h;
	m_f[0] = m_size[0];
	m_f[1] = m_size[1];
	m_c[0] = 0.5*m_f[0];
	m_c[1] = 0.5*m_f[1];

	for(size_t i=0;i<5;i++)
		m_k[i]=0;

	m_alpha = 0;

	m_F = mat(4,4);
	m_F.Eye();

	m_Finv = m_F;

}

ostream& operator<< (ostream& os, const CCam& x) {

	os << "# dims" << endl;
	os << x.m_size[0] << " " << x.m_size[1] << endl;
	os << "# focal lengths" << endl;
	os << x.m_f[0] << " " << x.m_f[1] << endl;
	os << "# principle point" << endl;
	os << x.m_c[0] << " " << x.m_c[1] << endl;
	os << "# radial distortion coefficients" << endl;
	os << x.m_k[0] << " " << x.m_k[1] << " " << x.m_k[2] << " " << x.m_k[3] << " " << x.m_k[4] << endl;
	os << "# skew coefficient" << endl;
	os << x.m_alpha << endl;
	os << "# frame world -> cam" << endl;

    for(size_t i=0; i<4; i++) {

        for(size_t j=0; j<4; j++) {
            os << x.m_F.Get(i,j);

            if(j<3)
                os << " ";
            else
                os << endl;

        }

    }

	return os;

}

bool CCam::SaveToFile(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file " << filename << "." << endl;
		return 1;

	 }

	out << *this;

	out.close();

	return 0;

}


bool CCam::OpenFromFile(const char* filename) {

	ifstream in(filename);

	if(!in) {

		cout << "ERROR: Could not open " << filename << "." << endl;
		return 1;

	 }


	string linebuffer;

	getline(in,linebuffer);

	in >> m_size[0];
	in >> m_size[1];
	in.get();

	getline(in,linebuffer);

	in >> m_f[0];
	in >> m_f[1];
	in.get();

	getline(in,linebuffer);

	in >> m_c[0];
	in >> m_c[1];
	in.get();

	getline(in,linebuffer);

	in >> m_k[0];
	in >> m_k[1];
	in >> m_k[2];
	in >> m_k[3];
	in >> m_k[4];
	in.get();

	getline(in,linebuffer);

	in >> m_alpha;
	in.get();

	getline(in,linebuffer);

    for(size_t i=0; i<3; i++) {

        for(size_t j=0; j<4; j++)
            in >> m_F(i,j);

        in.get();

    }

	in.close();

	m_Finv = CLinearAlgebra::InvertTransformation(m_F);

	return 0;

}


vec CCam::Project(vec x) {

	// result
	vec xp(2);

	// transform into camera coordinate system
	vec xc = CLinearAlgebra::TransformPoint(m_F,x);

	// check cheirality condition
	if(xc(2)<0) {

		xp(0) = -1000;
		xp(1) = -1000;
		return xp;

	}

	// pinhole projection
	vec xn(2);
	xn(0) = xc(0)/xc(2);
	xn(1) = xc(1)/xc(2);

	// radial distortion
	double r = xn.Norm2();

	vec dx(2);
	dx(0) = 2*m_k[2]*xn(0)*xn(1) + m_k[3]*(r*r + 2*xn(0)*xn(0));
	dx(1) = m_k[2]*(r*r + 2*xn(1)*xn(1)) + 2*m_k[3]*xn(0)*xn(1);

	double fac = 1 + m_k[0]*r*r + m_k[1]*r*r*r*r + m_k[4]*r*r*r*r*r*r;
	vec xd = xn*fac + dx;

	// transform to pixel coordinates
	xp(0) = m_f[0]*(xd(0) + m_alpha*xd(1)) + m_c[0];
	xp(1) = m_f[1]*xd(1) + m_c[1];

	return xp;

}


vec CCam::ProjectLocal(vec xc) {

	// pinhole projection
	vec xn(2);
	xn(0) = xc(0)/xc(2);
	xn(1) = xc(1)/xc(2);

	// radial distortion
	double r = xn.Norm2();

	vec dx(2);
	dx(0) = 2*m_k[2]*xn(0)*xn(1) + m_k[3]*(r*r + 2*xn(0)*xn(0));
	dx(1) = m_k[2]*(r*r + 2*xn(1)*xn(1)) + 2*m_k[3]*xn(0)*xn(1);

	double fac = 1 + m_k[0]*r*r + m_k[1]*r*r*r*r + m_k[4]*r*r*r*r*r*r;
	vec xd = xn*fac + dx;

	// transform to pixel coordinates
	vec xp(2);
	xp(0) = m_f[0]*(xd(0) + m_alpha*xd(1)) + m_c[0];
	xp(1) = m_f[1]*xd(1) + m_c[1];

	return xp;

}

vec CCam::ProjectDirectionLocal(vec dc) {

	// pinhole projection
	vec dn(2);
	dn(0) = dc(0)/dc(2);
	dn(1) = dc(1)/dc(2);

	// radial distortion
	double r = dn.Norm2();

	vec dd(2);
	dd(0) = 2*m_k[2]*dn(0)*dn(1) + m_k[3]*(r*r + 2*dn(0)*dn(0));
	dd(1) = m_k[2]*(r*r + 2*dn(1)*dn(1)) + 2*m_k[3]*dn(0)*dn(1);

	double fac = 1 + m_k[0]*r*r + m_k[1]*r*r*r*r + m_k[4]*r*r*r*r*r*r;

	vec dr = dn*fac + dd;

	// transform to pixel coordinates
	vec dp(2);
	dp(0) = m_f[0]*(dr(0) + m_alpha*dr(1));
	dp(1) = m_f[1]*dr(1);

	return dp;

}

void CCam::Project(vec x, vec& u, mat& J) {

	u = Project(x);

	J(0,0) = m_f[0]/x(2);
	J(0,2) = -(u(0)-m_c[0])/x(2);
	J(1,1) = m_f[1]/x(2);
	J(1,2) = -(u(1)-m_c[1])/x(2);

}

void CCam::ProjectLocal(vec xc, vec& u, mat& J) {

	u = ProjectLocal(xc);

	J(0,0) = m_f[0]/xc(2);
	J(0,2) = -(u(0)-m_c[0])/xc(2);
	J(1,1) = m_f[1]/xc(2);
	J(1,2) = -(u(1)-m_c[1])/xc(2);

}

vec CCam::Normalize(vec u) {

	// TODO: add undistortion

	vec xd(3);
	xd(0) = (1/m_f[0])*(u(0)-m_c[0]);
	xd(1) = (1/m_f[1])*(u(1)-m_c[1]);
	xd(2) = 1;

	return CLinearAlgebra::TransformDirectionBack(m_F,xd);

}

vec CCam::NormalizeLocal(const vec& u) const {

	vec xd(3);
	xd(0) = (1/m_f[0])*(u.Get(0)-m_c[0]);
	xd(1) = (1/m_f[1])*(u.Get(1)-m_c[1]);
	xd(2) = 1;

	return xd;

}

vec CCam::GetOrigin() {

	return CLinearAlgebra::GetTranslationVector(m_Finv);

}

mat CCam::GetProjectionMatrix() const {

	mat K(3,3);
	K(0,0) = m_f[0];
	K(0,1) = m_alpha;
	K(0,2) = m_c[0];
	K(1,1) = m_f[1];
	K(1,2) = m_c[1];
	K(2,2) = 1;

	return K;

}


mat CCam::GetInverseProjectionMatrix() const {

	mat Kinv(3,3);
	Kinv(0,0) = 1/m_f[0];
	Kinv(0,1) = -m_alpha/(m_f[0]*m_f[1]);
	Kinv(0,2) = (m_alpha*m_c[1])/(m_f[0]*m_f[1])-(m_c[0]/m_f[0]);
	Kinv(1,1) = 1/m_f[1];
	Kinv(1,2) = -(m_c[1]/m_f[1]);
	Kinv(2,2) = 1;

	return Kinv;

}


}
