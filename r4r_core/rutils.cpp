/*
 * utils.cpp
 *
 *  Created on: May 29, 2012
 *      Author: jbalzer
 */


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


mat CLinearAlgebra::Hat(vec omega) {

	mat Omega(3,3);

	Omega(0,1) = -omega.Get(2);
	Omega(0,2) = omega.Get(1);
	Omega(1,0) = omega.Get(2);
	Omega(1,2) = -omega.Get(0);
	Omega(2,0) = -omega.Get(1);
	Omega(2,1) = omega.Get(0);

	return Omega;

}

vec CLinearAlgebra::TransformPoint(mat& F, vec& x) {

	vec xh(x.NElems()+1);

	// homogenize
	for(size_t i=0;i<x.NElems();i++)
		xh(i) = x.Get(i);

	xh(x.NElems()) = 1;

	// transform
	xh = F*xh;

	vec result(x.NElems());

	// de-homogenize
	for(size_t i=0;i<x.NElems();i++)
		result(i) = xh(i);

	return result;

}


vec CLinearAlgebra::TransformDirection(mat& F, vec& x) {

	vec xh(x.NElems()+1);

	// homogenize
	for(size_t i=0;i<x.NElems();i++)
		xh(i) = x.Get(i);

	// transform
	xh = F*xh;

	vec result(x);

	// de-homogenize
	for(size_t i=0;i<x.NElems();i++)
		result(i) = xh(i);

	return result;

}

vec CLinearAlgebra::TransformPointBack(mat& F, vec& x) {

	mat Finv = InvertTransformation(F);

	return TransformPoint(Finv,x);

}

vec CLinearAlgebra::TransformDirectionBack(mat& F, vec& x) {

	mat Finv = InvertTransformation(F);

	return TransformDirection(Finv,x);

}

mat CLinearAlgebra::InvertTransformation(mat& F) {

	mat R(F.NRows()-1,F.NCols()-1);
	vec t(F.NRows()-1);

	for(size_t i=0; i<F.NRows()-1; i++) {

		for(size_t j=0; j<F.NCols()-1; j++) {

			R(i,j) = F(j,i);	// transpose!

		}

		t(i) = F(i,F.NCols()-1);

	}

	vec tinv = R*t;
	tinv.Scale(-1);

	// assemble result
	mat Finv(F);

	for(size_t i=0; i<Finv.NRows()-1; i++) {

		for(size_t j=0; j<Finv.NCols()-1; j++) {

			Finv(i,j) = R(i,j);

		}

		Finv(i,F.NCols()-1) = tinv(i);

	}

	return Finv;

}

mat CLinearAlgebra::Rodrigues(vec omega) {

	double theta = omega.Norm2();

	mat R(3,3);

    if(theta<std::numeric_limits<double>::epsilon()) {

		R(0,0) = 1;
		R(1,1) = 1;
		R(2,2) = 1;

		return R;

	}

	double ox, oy, oz;
	ox = omega(0)/theta;
	oy = omega(1)/theta;
	oz = omega(2)/theta;

	double oxox, oxoy, oxoz, oyoy, oyoz, ozoz;
	oxox = ox*ox;
	oxoy = ox*oy;
	oxoz = ox*oz;
	oyoy = oy*oy;
	oyoz = oy*oz;
	ozoz = oz*oz;

	double sth, cth, mcth;
    sth = sin(theta);
    cth  = cos(theta);
    mcth = 1 - cth;

    R(0,0) = 1 - mcth*(oyoy+ozoz);  R(0,1) = -sth*oz + mcth*oxoy;		R(0,2) = sth*oy + mcth*oxoz;
    R(1,0) = sth*oz + mcth*oxoy;	R(1,1) = 1 - mcth*(ozoz + oxox);	R(1,2) = -sth*ox + mcth*oyoz;
    R(2,0) = - sth*oy + mcth*oxoz;	R(2,1) = sth*ox + mcth*oyoz;    	R(2,2) = 1 - mcth*(oxox+oyoy);

    return R;

}

mat CLinearAlgebra::Rodrigues(double o1, double o2, double o3) {

	vec omega(3);
	omega(0) = o1;
	omega(1) = o2;
	omega(2) = o3;

	return Rodrigues(omega);

}

void CLinearAlgebra::Rodrigues(vec omega, mat& R, mat& DRx, mat& DRy, mat& DRz) {

	double theta = omega.Norm2();

    if(theta<std::numeric_limits<double>::epsilon()) {

	    R(0,0) = 1.0; 	R(0,1) = 0.0; 	R(0,2) = 0.0;
	    R(1,0) = 0.0; 	R(1,1) = 1.0; 	R(1,2) = 0.0;
	    R(2,0) = 0.0; 	R(2,1) = 0.0; 	R(2,2) = 1.0;

	    DRx(0,0) = 0; 	DRx(0,1) = 0;	DRx(0,2) = 0;
	    DRx(1,0) = 0;	DRx(1,1) = 0;   DRx(1,2) = -1;
	    DRx(2,0) = 0; 	DRx(2,1) = 1;	DRx(2,2) = 0;

	    DRy(0,0) = 0;	DRy(0,1) = 0;	DRy(0,2) = 1;
	    DRy(1,0) = 0;	DRy(1,1) = 0;	DRy(1,2) = 0;
	    DRy(2,0) = -1;	DRy(2,1) = 0;	DRy(2,2) = 0;

	    DRz(0,0) = 0;	DRz(0,1) = -1;	DRz(0,2) = 0;
	    DRz(1,0) = 1;	DRz(1,1) = 0;	DRz(1,2) = 0;
	    DRz(2,0) = 0;	DRz(2,1) = 0;	DRz(2,2) = 0;

	    return;

	}

	double ox, oy, oz;
	ox = omega(0)/theta;
	oy = omega(1)/theta;
	oz = omega(2)/theta;

	double oxox, oxoy, oxoz, oyoy, oyoz, ozoz;
	oxox = ox*ox;
	oxoy = ox*oy;
	oxoz = ox*oz;
	oyoy = oy*oy;
	oyoz = oy*oz;
	ozoz = oz*oz;

	double sth, cth, mcth;
    sth = sin(theta);
    cth  = cos(theta);
    mcth = 1 - cth;

    R(0,0) = 1 - mcth*(oyoy+ozoz);  R(0,1) = -sth*oz + mcth*oxoy;		R(0,2) = sth*oy + mcth*oxoz;
    R(1,0) = sth*oz + mcth*oxoy;	R(1,1) = 1 - mcth*(ozoz + oxox);	R(1,2) = -sth*ox + mcth*oyoz;
    R(2,0) = - sth*oy + mcth*oxoz;	R(2,1) = sth*ox + mcth*oyoz;    	R(2,2) = 1 - mcth*(oxox+oyoy);

    double a, b, c, d;
    a =  sth/theta;
    b = mcth/theta;
    c = cth - a;
    d = sth - 2*b;

    DRx(0,0) = -d*(oyoy + ozoz)*ox;			DRx(0,1) = b*oy - c*oxoz + d*oxoy*ox;		DRx(0,2) = b*oz + c*oxoy + d*oxoz*ox;
    DRx(1,0) = b*oy + c*oxoz + d*oxoy*ox;	DRx(1,1) = -2*b*ox - d*(ozoz + oxox)*ox;	DRx(1,2) = -a - c*oxox + d*oyoz*ox;
    DRx(2,0) = b*oz - c*oxoy + d*oxoz*ox; 	DRx(2,1) = a + c*oxox + d*oyoz*ox;			DRx(2,2) = -2*b*ox - d*(oyoy + oxox)*ox;

    DRy(0,0) = -2*b*oy - d*(oyoy + ozoz)*oy;	DRy(0,1) = b*ox - c*oyoz + d*oxoy*oy;	DRy(0,2) = a + c*oyoy + d*oxoz*oy;
    DRy(1,0) = b*ox + c*oyoz + d*oxoy*oy;		DRy(1,1) = -d*(ozoz + oxox)*oy;			DRy(1,2) =  b*oz - c*oxoy + d*oyoz*oy;
    DRy(2,0) = -a - c*oyoy + d*oxoz*oy;			DRy(2,1) =  b*oz + c*oxoy + d*oyoz*oy;	DRy(2,2) = -2*b*oy - d*(oyoy + oxox)*oy;

    DRz(0,0) = -2*b*oz - d*(oyoy + ozoz)*oz;	DRz(0,1) = -a - c*ozoz + d*oxoy*oz;			DRz(0,2) = b*ox + c*oyoz + d*oxoz*oz;
    DRz(1,0) = a + c*ozoz + d*oxoy*oz;			DRz(1,1) =  -2*b*oz - d*(ozoz + oxox)*oz;	DRz(1,2) = b*oy - c*oxoz + d*oyoz*oz;
    DRz(2,0) =  b*ox - c*oyoz + d*oxoz*oz;		DRz(2,1) = b*oy + c*oxoz + d*oyoz*oz;		DRz(2,2) = - d*(oyoy + oxox)*oz;

}

mat CLinearAlgebra::CreateTransformationMatrix(double tx, double ty, double tz, double ox, double oy, double oz) {

	vec omega(3);
	omega(0) = ox;
	omega(1) = oy;
	omega(2) = oz;

	mat R = Rodrigues(omega);

	mat F(4,4);
	F(3,3) = 1;

	F(0,3) = tx;
	F(1,3) = ty;
	F(2,3) = tz;

	for(size_t i=0; i<F.NRows()-1; i++) {

		for(size_t j=0; j<F.NCols()-1; j++) {

			F(i,j) = R(i,j);

		}

	}

	return F;

}

mat CLinearAlgebra::CreateInverseTransformationMatrix(double tx, double ty, double tz, double ox, double oy, double oz) {

	mat F = CreateTransformationMatrix(tx,ty,tz,ox,oy,oz);

	return InvertTransformation(F);

}


mat CLinearAlgebra::CreateTransformationMatrix(mat& R, vec& t) {

	mat F(4,4);
	F(3,3) = 1;

	for(size_t i=0; i<F.NRows()-1; i++) {

		for(size_t j=0; j<F.NCols()-1; j++) {

			F(i,j) = R(i,j);

		}

		F(i,3) = t(i);

	}

	return F;

}


mat CLinearAlgebra::GetRotationMatrix(mat& F) {

	mat R(F.NRows()-1,F.NCols()-1);

	for(size_t i=0; i<F.NRows()-1; i++) {

		for(size_t j=0; j<F.NCols()-1; j++) {

			R(i,j) = F(i,j);

		}

	}

	return R;

}

vec CLinearAlgebra::GetTranslationVector(mat& F) {

	vec t(F.NRows()-1);

	for(size_t i=0; i<F.NRows()-1; i++)
		t(i) = F(i,F.NCols()-1);

	return t;

}


double CLinearAlgebra::MedianAbsoluteDeviation(vec& x) {

	double median = x.Median();

	vec temp(x.NElems());

	for(size_t i=0; i<temp.NElems(); i++)
		temp(i) = fabs(x.Get(i)-median);

	return temp.Median();

}

mat CLinearAlgebra::ProjectToSO3(mat& R) {

	mat U(R.NRows(),R.NRows());
	vec s(min(R.NRows(),R.NCols()));
	mat Vt(R.NCols(),R.NCols());

	CMatrixFactorization::SVD(R,U,s,Vt);

	return U*Vt;

}

mat CLinearAlgebra::CalibratedNPoint(const vector<pair<vec,vec> >& corr, const CCam& cam) {

	mat E0(3,3);
	E0.Eye();

	if(corr.size()<8) {

		cout << "ERROR: Not enough point correspondences for estimation of essential matrix..." << endl;
		return E0;

	}

	mat A(corr.size(),9);

	for(size_t i=0; i<A.NRows(); i++) {

		vec x0 = cam.NormalizeLocal(corr[i].first);
		vec x1 = cam.NormalizeLocal(corr[i].second);

		A(i,0) = x0(0)*x1(0);
		A(i,1) = x0(0)*x1(1);
		A(i,2) = x0(0)*x1(2);
		A(i,3) = x0(1)*x1(0);
		A(i,4) = x0(1)*x1(1);
		A(i,5) = x0(1)*x1(2);
		A(i,6) = x0(2)*x1(0);
		A(i,7) = x0(2)*x1(1);
		A(i,8) = x0(2)*x1(2);

	}

	// solve least-squares problem by SVD
	mat U(A.NRows(),A.NRows());
	vec s(min(A.NRows(),A.NCols()));
	mat Vt(A.NCols(),A.NCols());

	// SVD
	CMatrixFactorization::SVD(A,U,s,Vt);

	// essential matrix builds from last row of Vt
	E0(0,0) = Vt(Vt.NRows()-1,0);
	E0(1,0) = Vt(Vt.NRows()-1,1);
    E0(2,0) = Vt(Vt.NRows()-1,2);
    E0(0,1) = Vt(Vt.NRows()-1,3);
    E0(1,1) = Vt(Vt.NRows()-1,4);
    E0(2,1) = Vt(Vt.NRows()-1,5);
    E0(0,2) = Vt(Vt.NRows()-1,6);
    E0(1,2) = Vt(Vt.NRows()-1,7);
    E0(2,2) = Vt(Vt.NRows()-1,8);

    // project to the space of essential matrices
	mat Up(3,3);
	vec sp(3);
	mat Vtp(3,3);

	// SVD
	CMatrixFactorization::SVD(E0,Up,sp,Vtp);

	double sm = 0.5*(sp(0)+sp(1));
	mat sigma(3,3);
	sigma(0,0) = sm;
	sigma(1,1) = sm;

	return Up*(sigma*Vtp);

}

mat CLinearAlgebra::EstimateHomography(const vector<pair<vec,vec> >& corr) {

	mat H(3,3);
	H.Eye();

	if(corr.size()<8) {

		cout << "ERROR: Not enough point correspondences for estimation of essential matrix..." << endl;
		return H;

	}


	mat A(2*corr.size(),9);

	for(size_t i=0; i<corr.size(); i++) {

		A(i,0) = corr[i].first.Get(0);
		A(i,1) = corr[i].first.Get(1);
		A(i,2) = 1;
		A(i,3) = 0;
		A(i,4) = 0;
		A(i,5) = 0;
		A(i,6) = -corr[i].second.Get(0)*corr[i].first.Get(0);
		A(i,7) = -corr[i].second.Get(0)*corr[i].first.Get(1);
		A(i,8) = -corr[i].second.Get(0);

		A(corr.size()+i,0) = 0;
		A(corr.size()+i,1) = 0;
		A(corr.size()+i,2) = 0;
		A(corr.size()+i,3) = corr[i].first.Get(0);
		A(corr.size()+i,4) = corr[i].first.Get(1);
		A(corr.size()+i,5) = 1;
		A(corr.size()+i,6) = -corr[i].second.Get(1)*corr[i].first.Get(0);
		A(corr.size()+i,7) = -corr[i].second.Get(1)*corr[i].first.Get(1);
		A(corr.size()+i,8) = -corr[i].second.Get(1);

	}

	mat U(A.NRows(),A.NRows());
	vec s(min(A.NRows(),A.NCols()));
	mat Vt(A.NCols(),A.NCols());

	// SVD
	CMatrixFactorization::SVD(A,U,s,Vt);

	// assemle homography
	H(0,0) = Vt(Vt.NRows()-1,0);
	H(0,1) = Vt(Vt.NRows()-1,1);
    H(0,2) = Vt(Vt.NRows()-1,2);
    H(1,0) = Vt(Vt.NRows()-1,3);
    H(1,1) = Vt(Vt.NRows()-1,4);
    H(1,2) = Vt(Vt.NRows()-1,5);
    H(2,0) = Vt(Vt.NRows()-1,6);
    H(2,1) = Vt(Vt.NRows()-1,7);
    H(2,2) = Vt(Vt.NRows()-1,8);

  	// normalize
    H.Scale(1/H.Get(2,2));
    //if(H(2,2)<0)
    //	H.Scale(-1);

	return H;

}



mat CLinearAlgebra::FactorEssentialMatrix(const mat& E) {

	mat En(E);
	En.Normalize();

	mat Rzp(3,3);
	Rzp(0,1) = 1;
	Rzp(1,0) = -1;
	Rzp(2,2) = 1;

	mat Rzm = mat::Transpose(Rzp);

	mat U(3,3);
	vec s(3);
	mat Vt(3,3);
	mat Ut;
	mat S(3,3);

	// init result
	mat R;
	vec t(3);

	// first try with positive sign, there are two options
	CMatrixFactorization::SVD(En,U,s,Vt);

	Ut = mat::Transpose(U);

	for(size_t i=0; i<3; i++)
		S(i,i)=s(i);

	R = U*(Rzm*Vt);

	if(R.Determinant()>0 && R(2,2)>0) {

		mat T = U*(Rzp*(S*Ut));

		t(0) = -T(1,2);
		t(1) = T(0,2);
		t(2) = -T(0,1);

		return CreateTransformationMatrix(R,t);

	}

	// flipped configuration
	R = U*(Rzp*Vt);

	if(R.Determinant()>0 && R(2,2)>0) {

		mat T = U*(Rzm*(S*Ut));

		t(0) = -T(1,2);
		t(1) = T(0,2);
		t(2) = -T(0,1);

		return CreateTransformationMatrix(R,t);

	}

	// try with different sign
	En.Scale(-1);

	CMatrixFactorization::SVD(En,U,s,Vt);

	Ut = mat::Transpose(U);

	for(size_t i=0; i<3; i++)
		S(i,i)=s(i);

	R = U*(Rzm*Vt);

	if(R.Determinant()>0 && R(2,2)>0) {

		mat T = U*(Rzp*(S*Ut));

		t(0) = -T(1,2);
		t(1) = T(0,2);
		t(2) = -T(0,1);

		return CreateTransformationMatrix(R,t);

	}

	R = U*(Rzp*Vt);

	if(R.Determinant()>0 && R(2,2)>0) {

		mat T = U*(Rzm*(S*Ut));

		t(0) = -T(1,2);
		t(1) = T(0,2);
		t(2) = -T(0,1);

		return CreateTransformationMatrix(R,t);

	}

	// otherwise return identity transformation
	R.Eye();
	t.Scale(0);

	cout << "WARNING: Factorization failed due to ambiguity..." << endl;

	return CreateTransformationMatrix(R,t);

}

mat CLinearAlgebra::ZhangFactorization(const mat& H, const CCam& cam) {

	mat Kinv = cam.GetInverseProjectionMatrix();

	// rotational part
	vec ex =  Kinv*(H.GetColumn(0));
	double l = ex.Norm2();
	ex.Scale(1.0/l);
	vec ey = Kinv*(H.GetColumn(1));
	ey.Normalize();
	vec ez = vec::CrossProduct(ex,ey);

	mat R(3,3);
	R.SetColumn(0,ex);
	R.SetColumn(1,ey);
	R.SetColumn(2,ez);

	// translational part
	vec t = Kinv*(H.GetColumn(2));
	t.Scale(1.0/l);

	// project onto SO(3)
	R = CLinearAlgebra::ProjectToSO3(R);

	// assemble frame
	mat F = CLinearAlgebra::CreateTransformationMatrix(R,t);

	return F;

}





}
