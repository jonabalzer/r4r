/*
 * trafo.cpp
 *
 *  Created on: Oct 16, 2012
 *      Author: jbalzer
 */




#include "trafo.h"
#include "rutils.h"
#include "factor.h"
#include <math.h>
#include <assert.h>

using namespace std;

namespace R4R {

template <size_t size>
CTransformation<size>::CTransformation():
	m_F(size+1,size+1) {

	m_F.Eye();

}

template <size_t size>
CTransformation<size>::CTransformation(const mat& F) {

	assert(F.NRows()==size+1 && F.NCols()==size+1);

	m_F = mat(F);

}

template <size_t size>
mat CTransformation<size>::GetLinearPart() const {

	mat A(size,size);

	for(size_t i=0; i<size; i++) {

		for(size_t j=0; j<size; j++) {

			A(i,j) = m_F.Get(i,j);

		}

	}

	return A;

}


template <size_t size>
vec CTransformation<size>::GetTranslation() const {

	vec t(size);

	for(size_t i=0; i<size; i++)
		t(i) = m_F.Get(i,size);

	return t;

}

template <size_t size>
void CTransformation<size>::SetLinearPart(const mat& A) {

	for(size_t i=0; i<size; i++) {

		for(size_t j=0; j<size; j++) {

			m_F(i,j) = A.Get(i,j);

		}

	}

}

template <size_t size>
void CTransformation<size>::SetTranslation(const vec& t) {

	for(size_t i=0; i<size; i++)
		m_F(i,size) = t.Get(i);

}



template class CTransformation<2>;
template class CTransformation<3>;


template <size_t size>
CRotation<size>::CRotation(const mat& A) {

	assert(A.NRows()==size && A.NCols()==size);

	mat U(A.NRows(),A.NRows());
	vec s(min(A.NRows(),A.NCols()));
	mat Vt(A.NCols(),A.NCols());

    CMatrixFactorization<double>::SVD(A,U,s,Vt);

	mat Ut = mat::Transpose(U);

	mat R = U*Ut;

}

mat CRotation<2>::Rodrigues(double o) {

	mat R(2,2);

	R(0,0) = cos(o);
	R(0,1) = -sin(o);
	R(1,0) = sin(o);
	R(1,1) = cos(o);

	return R;

}




mat CRotation<3>::Rodrigues(double o1, double o2, double o3) {

	double theta = sqrt(o1*o1+o2*o2+o3*o3);

	mat R(3,3);

	if(theta<1e-20) {

		R(0,0) = 1;
		R(1,1) = 1;
		R(2,2) = 1;

		return R;

	}

	double ox, oy, oz;
	ox = o1/theta;
	oy = o2/theta;
	oz = o3/theta;

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

CRotation<3>::CRotation(double o1, double o2, double o3) {

	mat R = CRotation<3>::Rodrigues(o1,o2,o3);

	CTransformation<3>::SetLinearPart(R);

}

CRotation<3>::CRotation(vec omega) {

	assert(omega.NElems()==3);

	CRotation(omega.Get(0),omega.Get(1),omega.Get(2));

}

vec CRotation<3>::Log(const mat& R) {

	vec omega(3);

	double arg = 0.5*(R.Trace()-1.0);

	if(fabs(arg)>=1)
		return omega;

	double theta = acos(arg);

	mat Omega = R - mat::Transpose(R);
	Omega.Scale(theta/(2*sin(theta)));

	omega(0) = -Omega.Get(1,2);
	omega(1) = Omega.Get(0,2);
	omega(2) = -Omega.Get(0,1);

	return omega;

}

vec CRotation<3>::Log() {

	mat R = GetLinearPart();

	return Log(R);

}


template class CRotation<2>;
template class CRotation<3>;



}
