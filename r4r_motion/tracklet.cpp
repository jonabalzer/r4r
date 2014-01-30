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
#include <fstream>
#include <limits>
#include <algorithm>

#include "descriptor.h"
#include "tracklet.h"

using namespace cv;
using namespace std;

namespace R4R {

template<template<class T, class Allocator = std::allocator<T> > class Container>
CTracklet<Container>::CTracklet(size_t t0, const imfeature& x0, size_t maxlength):
    m_data(maxlength),
	m_t0(t0),
    m_status(true),
    m_hash(CTracklet::GenerateHash(t0,x0)) {

    m_data.push_back(x0);

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracklet<Container>::Update(const imfeature& x) {

    m_data.push_back(x);

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
const vec2f& CTracklet<Container>::GetPastLocation(size_t steps) const {

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
vec2f CTracklet<Container>::GetPastLocationAtNativeScale(size_t steps) const {

}

template<template<class T, class Allocator = std::allocator<T> > class C>
ostream& operator<<(ostream& os, const CTracklet<C>& x) {


 //   os << x.m_data.at(i%x.m_data.size()) << endl;

    //return os;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
string CTracklet<Container>::GenerateHash(size_t t0, const imfeature& x) {

    stringstream t;
    t.fill('0');
    t.width(4);
    t << t0;

    const vec2f& pt = x.GetLocation();

    stringstream x0, y0;
    x0.fill('0');
    x0.width(4);
    x0 << size_t(pt.Get(0)+0.5);
    y0.fill('0');
    y0.width(4);
    y0 << size_t(pt.Get(1)+0.5);

    stringstream scale;
    scale.fill('0');
    scale.width(4);
    scale << size_t(x.GetScale());

    stringstream hash;
    hash << t.str() << x0.str() << y0.str() << scale.str();

    return hash.str();

}

#ifdef QT_GUI_LIB
template<template<class T, class Allocator = std::allocator<T> > class Container>
const Qt::GlobalColor CTracklet<Container>::COLORS[10] = { Qt::green,
                                                           Qt::red,
                                                           Qt::blue,
                                                           Qt::cyan,
                                                           Qt::magenta,
                                                           Qt::white,
                                                           Qt::yellow,
                                                           Qt::darkGreen,
                                                           Qt::darkRed,
                                                           Qt::darkBlue };

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracklet<Container>::Draw(QImage& img, size_t length) const {

    // check if length > 0
 /*   if(length==0)
        return;

    // get latest feature
    imfeature ft = m_data.at(m_cursor);

    // get scale for color, mod 10 because so far only 10 colors are available
    u_int scale = u_int(ft.GetScale())%10;

    // prepare drawing
    QPainter p(&img);
    QPen pen;
    pen.setColor(COLORS[scale]);
    pen.setWidth(2);
    p.setPen(pen);

    // draw circle if length of trail is 1
    if(length==1) {

        vec2f xt = ft.GetLocationAtNativeScale();
        p.drawEllipse(int(xt.Get(0)),int(xt.Get(1)),10,10);
        return;

    }

    // draw trail if desired, not longer than length or size of container
    for(size_t i=1; i<length && i<m_data.size(); i++) {

        // get next feature
        const imfeature& ftm1 = m_data.at((m_cursor+i)%m_data.size());

        // extract locations
        vec2f xt = ft.GetLocationAtNativeScale();
        vec2f xtm1 = ftm1.GetLocationAtNativeScale();

        // only draw if the features do not come from the tracklet initialization
        if(!xtm1.IsZero() && !xt.IsZero())
            p.drawLine(int(xt.Get(0)),int(xt.Get(1)),int(xtm1.Get(0)),int(xtm1.Get(1)));
        else
            break;

        // make the current feature, the "oldest"
        ft = ftm1;

    }*/

}
#endif

template class CTracklet<list>;
template class CTracklet<vector>;


} // end of namespace



