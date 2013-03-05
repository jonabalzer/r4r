/*
 * rect.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: jbalzer
 */

#include "rect.h"
#include <iostream>
#include <assert.h>

#define MARGIN 2


using namespace cv;
using namespace std;


namespace R4R {

template <class T>
CRectangle<T>::CRectangle() {

	m_t[0] = 0;
	m_t[1] = 0;
	m_s[0] = 1;
	m_s[1] = 1;
	m_phi = 0;

}


template <class T>
CRectangle<T>::CRectangle(T tx, T ty, T sx, T sy, double phi) {

	m_t[0] = tx;
	m_t[1] = ty;
	m_s[0] = sx;
	m_s[1] = sy;
	m_phi = phi;

}

template <class T>
CDenseVector<T> CRectangle<T>::TransformFrom(T x, T y) {

	CDenseVector<T> result(2);

	if(m_phi!=0) {

		T cphi = (T)cos(m_phi);
		T sphi = (T)sin(m_phi);

		result(0) = m_s[0]*(cphi*x - sphi*y) + m_t[0];
		result(1) = m_s[1]*(sphi*x + cphi*y) + m_t[1];

	} else {

		result(0) = m_s[0]*x + m_t[0];
		result(1) = m_s[1]*y + m_t[1];

	}

	return result;

}


template <class T>
CDenseVector<T> CRectangle<T>::Barycenter() {

	CDenseVector<T> result(2);

	result(0) = m_t[0];
	result(1) = m_t[1];

	return result;

}

template <class T>
void CRectangle<T>::Draw(cv::Mat& img, cv::Scalar color) {

	if(m_phi!=0) {

		CDenseVector<T> tl = TopLeft();
		CDenseVector<T> bl = BottomLeft();
		CDenseVector<T> br = BottomRight();
		CDenseVector<T> tr = TopRight();

		line(img,Point2f(tl(0),tl(1)),Point2f(bl(0),bl(1)),color);
		line(img,Point2f(bl(0),bl(1)),Point2f(br(0),br(1)),color);
		line(img,Point2f(br(0),br(1)),Point2f(tr(0),tr(1)),color);
		line(img,Point2f(tr(0),tr(1)),Point2f(tl(0),tl(1)),color);

	}
	else {

		CDenseVector<T> tl = TopLeft();
		CDenseVector<T> br = BottomRight();

		rectangle(img,Point2f(tl(0),tl(1)),Point2f(br(0),br(1)),color);


	}


}


template <class U>
ostream& operator<<(ostream& os, CRectangle<U>& x) {

	vec tl = x.TopLeft();
	vec bl = x.BottomLeft();
	vec br = x.BottomRight();
	vec tr = x.TopRight();

	tl.Transpose();
	bl.Transpose();
	br.Transpose();
	tr.Transpose();

	os << tl << endl;
	os << bl << endl;
	os << br << endl;
	os << tr << endl;

	return os;

}


template class CRectangle<size_t>;
template class CRectangle<float>;
template class CRectangle<double>;
template ostream& operator<< (ostream& os, CRectangle<double>& x);

}
