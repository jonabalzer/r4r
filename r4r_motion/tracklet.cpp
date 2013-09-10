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

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "descriptor.h"
#include "tracklet.h"

using namespace cv;
using namespace std;

namespace R4R {

CTracklet::CTracklet():
    list<imfeature>(),
	m_t0(0),
	m_status(false) {}

CTracklet::CTracklet(size_t t0, imfeature x0):
    list<imfeature>(),
	m_t0(t0),
	m_status(true)
{

	push_back(x0);

}

void CTracklet::Update(imfeature x) {

	push_back(x);

}

vec2f CTracklet::GetLatestVelocity() {

    list<imfeature>::reverse_iterator rit = rbegin();

    vec2f x1 = (*rit).GetLocation();

    vec2f v;

	if(size()>1) {

        rit++;

        vec2f x0 = (*rit).GetLocation();

        v = x1 - x0;

	}

	return v;

}

#ifdef QT_GUI_LIB

const Qt::GlobalColor CTracklet::COLORS[10] = { Qt::green,
                                                Qt::red,
                                                Qt::blue,
                                                Qt::cyan,
                                                Qt::magenta,
                                                Qt::white,
                                                Qt::yellow,
                                                Qt::darkGreen,
                                                Qt::darkRed,
                                                Qt::darkBlue };

void CTracklet::Draw(QImage& img, size_t length) const {

    // check if length > 0
    if(length==0)
        return;

    // access to the two last features
    list<imfeature>::const_reverse_iterator ita, itb;
    ita = rbegin();
    itb = rbegin();
    itb++;

    // get scale for color
    u_int scale = u_int((*ita).GetScale());

    // prepare drawing
    QPainter p(&img);
    QPen pen;
    pen.setColor(COLORS[scale]);
    pen.setWidth(2);
    p.setPen(pen);

    // draw circle if length of trail is 1
    if(length==1) {

        vec2f a = (*ita).GetLocationAtNativeScale();
        p.drawEllipse(int(a.Get(0)),int(a.Get(1)),10,10);
        return;

    }

    size_t counter = 0;

    while (itb!=rend() && counter<length) {

        vec2f a = (*ita).GetLocationAtNativeScale();
        vec2f b = (*itb).GetLocationAtNativeScale();

        p.drawLine(int(a.Get(0)),int(a.Get(1)),int(b.Get(0)),int(b.Get(1)));

    	ita++;
    	itb++;
    	counter++;

    }

}

#endif

ostream& operator<<(ostream& os, CTracklet& x) {

    list<imfeature>::iterator it;

	for(it=x.begin(); it!=x.end(); it++)
        os << (*it) << endl;

	return os;

}

std::string CTracklet::GetHash() const {

	stringstream t0;
	t0.fill('0');
	t0.width(4);
	t0 << m_t0;

    vec2f pt = front().GetLocation();

	stringstream x0, y0;
	x0.fill('0');
	x0.width(4);
    x0 << (size_t)(pt.Get(0)+0.5);
	y0.fill('0');
	y0.width(4);
    y0 << (size_t)(pt.Get(1)+0.5);

	stringstream hash;
	hash << t0.str() << x0.str() << y0.str();

	return hash.str();

}

vec2f CTracklet::GetPastLocation(size_t steps) const {

    list<imfeature>::const_reverse_iterator rit = rbegin();

	for(size_t i=0; i<steps; i++)
		rit++;

    return rit->GetLocation();

}

vec2f CTracklet::GetPastLocationAtNativeScale(size_t steps) const {

    list<imfeature>::const_reverse_iterator rit = rbegin();

	for(size_t i=0; i<steps; i++)
		rit++;

    return rit->GetLocationAtNativeScale();

}

}



