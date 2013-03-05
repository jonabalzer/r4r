/*
 * tracklet.cpp
 *
 *  Created on: Mar 16, 2012
 *      Author: jbalzer
 */

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "descriptor.h"
#include "tracklet.h"

using namespace cv;
using namespace std;

namespace R4R {

CTracklet::CTracklet():
    list<CFeature>(),
	m_t0(0),
	m_scale(0),
	m_status(false) {

}

CTracklet::CTracklet(size_t t0, size_t s, CFeature x0):
    list<CFeature>(),
	m_t0(t0),
	m_scale(s),
	m_status(true)
{

	push_back(x0);

}

CTracklet::~CTracklet() {

//	list<shared_ptr<CFeature> >::iterator it;

//	for(it=begin(); it!=end(); it++)
//		it->reset();

	clear();

}


void CTracklet::Update(CFeature x) {

	push_back(x);

}


vec CTracklet::GetLatestVelocity() {

    list<CFeature>::reverse_iterator rit = rbegin();

    vec x1 = (*rit).GetLocation();

	vec v(2);

	if(size()>1) {

	rit++;

    vec x0 = (*rit).GetLocation();

	v = x1 - x0;

	}

	return v;

}


void CTracklet::Draw(Mat& img, size_t length) {

    list<CFeature>::reverse_iterator ita, itb;

    ita = rbegin();
    itb = rbegin();
    itb++;

    size_t counter = 0;

    while (itb!=rend() && counter<length) {

        vec a = (*ita).GetLocation();
        vec b = (*itb).GetLocation();

    	double s = pow(2,m_scale);
    	a.Scale(s);
    	b.Scale(s);

    	line(img,Point2f(a(0),a(1)),Point2f(b(0),b(1)),CFeature::COLORS[m_scale]);

    	ita++;
    	itb++;
    	counter++;

    }

}

ostream& operator<<(ostream& os, CTracklet& x) {

    list<CFeature>::iterator it;

	for(it=x.begin(); it!=x.end(); it++)
        os << (*it);

	return os;

}

bool CTracklet::SaveToFile(const char* filename) {

    CFeature::SaveToFile(filename,*this);

//	if(size()==0)
//		return 1;

//	ofstream out(filename,ofstream::binary);

//	if(!out) {

//		cout << "ERROR: Could not open file " << filename << "." << endl;
//		return 1;

//	 }

//	out << "# Creation time:" << endl << m_t0 << endl;
//	out << "# Life time:" << endl << size() << endl;
//	out << "# Scale:" << endl << m_scale << endl;
//	out << "# x-coordinate y-coordinate" << endl;
//	out << "# quality" << endl;
//	out << "# number of descriptor" << endl;
//	out <<	"# descriptor 1" << endl;
//	out <<	"# ..." << endl;
//	out <<	"# descriptor n" << endl;

//	out << *this;

//	out.close();

	return 0;
}

std::string CTracklet::GetHash() {

	stringstream t0;
	t0.fill('0');
	t0.width(4);
	t0 << m_t0;

    vec pt = front().GetLocation();

	stringstream x0, y0;
	x0.fill('0');
	x0.width(4);
	x0 << pt.Get(0);
	y0.fill('0');
	y0.width(4);
	y0 << pt.Get(1);

	stringstream hash;
	hash << t0.str() << x0.str() << y0.str();

	return hash.str();

}

void CTracklet::DeleteDescriptors() {

    list<CFeature>::iterator it;

	for(it=begin(); it!=end(); it++)
        it->DeleteDescriptors();

}

vec CTracklet::ComputeVariance() {

	vec x(size()), y(size());
	vec result(2);

    list<CFeature>::iterator it;

	size_t t = 0;

	for(it=begin(); it!=end(); it++) {

        vec location = it->GetLocation();

		x(t) = location(0);
		y(t) = location(1);

		t++;

	}

	result(0) = x.Variance();
	result(1) = y.Variance();

	return result;

}

vec CTracklet::GetPastLocation(size_t steps) {

    list<CFeature>::reverse_iterator rit = rbegin();

	for(size_t i=0; i<steps; i++)
		rit++;

    return rit->GetLocation();

}

vec CTracklet::GetPastLocationAtNativeScale(size_t steps) {

    list<CFeature>::reverse_iterator rit = rbegin();

	for(size_t i=0; i<steps; i++)
		rit++;

    return rit->GetLocationAtNativeScale();

}

}



