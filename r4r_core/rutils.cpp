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

#include "rutils.h"
#include "factor.h"

#include <math.h>
#include <assert.h>
#include <limits>

using namespace std;

namespace R4R {

template<class T>
void CLinearAlgebra::GivensRotation(T a, T b, T& c, T& s, T& r) {

	if(b==0) {

		c = copysign(1.0,a);
		s = 0;
		r = fabs(a);

	}

	if(a==0) {

		c = 0;
		s = copysign(1.0,b);
		r = fabs(b);

	}

	if(fabs(b)>fabs(a)) {

		double t = a/b;
		double u = copysign(sqrt(1+t*t),b);
		s = 1/u;
		c = s*t;
		r = b*u;

	}
	else {

		double t = b/a;
		double u = copysign(sqrt(1.0+t*t),a);
		c = 1/u;
		s = c*t;
		r = a*u;

	}

}

template void CLinearAlgebra::GivensRotation(double a, double b, double& c, double& s, double& r);
template void CLinearAlgebra::GivensRotation(float a, float b, float& c, float& s, float& r);

mat CLinearAlgebra::CalibratedNPoint(const vector<pair<vec,vec> >& corr, const CPinholeCam<double>& cam) {

//	mat E0(3,3);
//	E0.Eye();

//	if(corr.size()<8) {

//		cout << "ERROR: Not enough point correspondences for estimation of essential matrix..." << endl;
//		return E0;

//	}

//	mat A(corr.size(),9);

//	for(size_t i=0; i<A.NRows(); i++) {

//		vec x0 = cam.NormalizeLocal(corr[i].first);
//		vec x1 = cam.NormalizeLocal(corr[i].second);

//		A(i,0) = x0(0)*x1(0);
//		A(i,1) = x0(0)*x1(1);
//		A(i,2) = x0(0)*x1(2);
//		A(i,3) = x0(1)*x1(0);
//		A(i,4) = x0(1)*x1(1);
//		A(i,5) = x0(1)*x1(2);
//		A(i,6) = x0(2)*x1(0);
//		A(i,7) = x0(2)*x1(1);
//		A(i,8) = x0(2)*x1(2);

//	}

//	// solve least-squares problem by SVD
//	mat U(A.NRows(),A.NRows());
//	vec s(min(A.NRows(),A.NCols()));
//	mat Vt(A.NCols(),A.NCols());

//	// SVD
//    CMatrixFactorization<double>::SVD(A,U,s,Vt);

//	// essential matrix builds from last row of Vt
//	E0(0,0) = Vt(Vt.NRows()-1,0);
//	E0(1,0) = Vt(Vt.NRows()-1,1);
//    E0(2,0) = Vt(Vt.NRows()-1,2);
//    E0(0,1) = Vt(Vt.NRows()-1,3);
//    E0(1,1) = Vt(Vt.NRows()-1,4);
//    E0(2,1) = Vt(Vt.NRows()-1,5);
//    E0(0,2) = Vt(Vt.NRows()-1,6);
//    E0(1,2) = Vt(Vt.NRows()-1,7);
//    E0(2,2) = Vt(Vt.NRows()-1,8);

//    // project to the space of essential matrices
//	mat Up(3,3);
//	vec sp(3);
//	mat Vtp(3,3);

//	// SVD
//    CMatrixFactorization<double>::SVD(E0,Up,sp,Vtp);

//	double sm = 0.5*(sp(0)+sp(1));
//	mat sigma(3,3);
//	sigma(0,0) = sm;
//	sigma(1,1) = sm;

//	return Up*(sigma*Vtp);

}

mat CLinearAlgebra::EstimateHomography(const vector<pair<vec,vec> >& corr) {

//    mat H(3,3);
//    H.Eye();

//    if(corr.size()<8) {

//        cout << "ERROR: Not enough point correspondences for estimation of essential matrix..." << endl;
//        return H;

//    }


//    mat A(2*corr.size(),9);

//    for(size_t i=0; i<corr.size(); i++) {

//        A(i,0) = corr[i].first.Get(0);
//        A(i,1) = corr[i].first.Get(1);
//        A(i,2) = 1;
//        A(i,3) = 0;
//        A(i,4) = 0;
//        A(i,5) = 0;
//        A(i,6) = -corr[i].second.Get(0)*corr[i].first.Get(0);
//        A(i,7) = -corr[i].second.Get(0)*corr[i].first.Get(1);
//        A(i,8) = -corr[i].second.Get(0);

//        A(corr.size()+i,0) = 0;
//        A(corr.size()+i,1) = 0;
//        A(corr.size()+i,2) = 0;
//        A(corr.size()+i,3) = corr[i].first.Get(0);
//        A(corr.size()+i,4) = corr[i].first.Get(1);
//        A(corr.size()+i,5) = 1;
//        A(corr.size()+i,6) = -corr[i].second.Get(1)*corr[i].first.Get(0);
//        A(corr.size()+i,7) = -corr[i].second.Get(1)*corr[i].first.Get(1);
//        A(corr.size()+i,8) = -corr[i].second.Get(1);

//    }

//    mat U(A.NRows(),A.NRows());
//    vec s(min(A.NRows(),A.NCols()));
//    mat Vt(A.NCols(),A.NCols());

//    // SVD
//    CMatrixFactorization<double>::SVD(A,U,s,Vt);

//    // assemle homography
//    H(0,0) = Vt(Vt.NRows()-1,0);
//    H(0,1) = Vt(Vt.NRows()-1,1);
//    H(0,2) = Vt(Vt.NRows()-1,2);
//    H(1,0) = Vt(Vt.NRows()-1,3);
//    H(1,1) = Vt(Vt.NRows()-1,4);
//    H(1,2) = Vt(Vt.NRows()-1,5);
//    H(2,0) = Vt(Vt.NRows()-1,6);
//    H(2,1) = Vt(Vt.NRows()-1,7);
//    H(2,2) = Vt(Vt.NRows()-1,8);

//    // normalize
//    H.Scale(1/H.Get(2,2));
//    //if(H(2,2)<0)
//    //	H.Scale(-1);
//
//    return H;

}



mat CLinearAlgebra::FactorEssentialMatrix(const mat& E) {

//	mat En(E);
//	En.Normalize();

//	mat Rzp(3,3);
//	Rzp(0,1) = 1;
//	Rzp(1,0) = -1;
//	Rzp(2,2) = 1;

//	mat Rzm = mat::Transpose(Rzp);

//	mat U(3,3);
//	vec s(3);
//	mat Vt(3,3);
//	mat Ut;
//	mat S(3,3);

//	// init result
//	mat R;
//	vec t(3);

//	// first try with positive sign, there are two options
//    CMatrixFactorization<double>::SVD(En,U,s,Vt);

//	Ut = mat::Transpose(U);

//	for(size_t i=0; i<3; i++)
//		S(i,i)=s(i);

//	R = U*(Rzm*Vt);

//	if(R.Determinant()>0 && R(2,2)>0) {

//		mat T = U*(Rzp*(S*Ut));

//		t(0) = -T(1,2);
//		t(1) = T(0,2);
//		t(2) = -T(0,1);

//		return CreateTransformationMatrix(R,t);

//	}

//	// flipped configuration
//	R = U*(Rzp*Vt);

//	if(R.Determinant()>0 && R(2,2)>0) {

//		mat T = U*(Rzm*(S*Ut));

//		t(0) = -T(1,2);
//		t(1) = T(0,2);
//		t(2) = -T(0,1);

//		return CreateTransformationMatrix(R,t);

//	}

//	// try with different sign
//	En.Scale(-1);

//    CMatrixFactorization<double>::SVD(En,U,s,Vt);

//	Ut = mat::Transpose(U);

//	for(size_t i=0; i<3; i++)
//		S(i,i)=s(i);

//	R = U*(Rzm*Vt);

//	if(R.Determinant()>0 && R(2,2)>0) {

//		mat T = U*(Rzp*(S*Ut));

//		t(0) = -T(1,2);
//		t(1) = T(0,2);
//		t(2) = -T(0,1);

//		return CreateTransformationMatrix(R,t);

//	}

//	R = U*(Rzp*Vt);

//	if(R.Determinant()>0 && R(2,2)>0) {

//		mat T = U*(Rzm*(S*Ut));

//		t(0) = -T(1,2);
//		t(1) = T(0,2);
//		t(2) = -T(0,1);

//		return CreateTransformationMatrix(R,t);

//	}

//	// otherwise return identity transformation
//	R.Eye();
//	t.Scale(0);

//	cout << "WARNING: Factorization failed due to ambiguity..." << endl;

//	return CreateTransformationMatrix(R,t);

}

mat CLinearAlgebra::ZhangFactorization(const mat& H, const CPinholeCam<double> &cam) {

//	mat Kinv = cam.GetInverseProjectionMatrix();

//	// rotational part
//	vec ex =  Kinv*(H.GetColumn(0));
//	double l = ex.Norm2();
//	ex.Scale(1.0/l);
//	vec ey = Kinv*(H.GetColumn(1));
//	ey.Normalize();
//	vec ez = vec::CrossProduct(ex,ey);

//	mat R(3,3);
//	R.SetColumn(0,ex);
//	R.SetColumn(1,ey);
//	R.SetColumn(2,ez);

//	// translational part
//	vec t = Kinv*(H.GetColumn(2));
//	t.Scale(1.0/l);

//	// project onto SO(3)
//	R = CLinearAlgebra::ProjectToSO3(R);

//	// assemble frame
//	mat F = CLinearAlgebra::CreateTransformationMatrix(R,t);

//	return F;

}





}
