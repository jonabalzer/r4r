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
    m_data(),
	m_t0(t0),
    m_status(true),
    m_hash(GenerateHash(t0,x0)) {

    m_data.push_back(x0);

}

template<>
CTracklet<CRingBuffer>::CTracklet(size_t t0, const imfeature& x0, size_t maxlength):
    m_data(maxlength),
    m_t0(t0),
    m_status(true),
    m_hash(GenerateHash(t0,x0)) {

    m_data.push_back(x0);

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracklet<Container>::Update(const imfeature& x) {

    m_data.push_back(x);

}

template<template<class T,class Allocator = std::allocator<T> > class Container>
const vec2f& CTracklet<Container>::GetPastLocation(size_t steps) const {

    typename Container<imfeature>::const_reverse_iterator rit = m_data.rbegin();

    for(size_t i=0; i<steps; i++)
        ++rit;

    return rit->GetLocation();

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
vec2f CTracklet<Container>::GetPastLocationAtNativeScale(size_t steps) const {

    typename Container<imfeature>::const_reverse_iterator rit = m_data.rbegin();

    for(size_t i=0; i<steps; i++)
        ++rit;

    return rit->GetLocationAtNativeScale();

}

template<template<class T, class Allocator = std::allocator<T> > class C>
ostream& operator<<(ostream& os, const CTracklet<C>& x) {

    os << "Hash: " << x.m_hash << endl;
    os << x.m_data.size() << endl;

    typename C<imfeature>::const_iterator it;

    for(it=x.m_data.begin(); it!=x.m_data.end(); ++it)
        os << (*it) << endl;

    return os;

}

template ostream& operator<<(ostream& os, const CTracklet<CRingBuffer>& x);
template ostream& operator<<(ostream& os, const CTracklet<list>& x);
template ostream& operator<<(ostream& os, const CTracklet<vector>& x);

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

    typename Container<imfeature>::const_reverse_iterator ita = m_data.rbegin();
    u_int scale = u_int((*ita).GetScale());

    // prepare drawing
    QPainter p(&img);
    QPen pen;
    pen.setColor(COLORS[scale]);
    pen.setWidth(2);
    p.setPen(pen);

    // draw circle if length of trail is 1
    if(length==0) {

        vec2f a = (*ita).GetLocationAtNativeScale();
        p.drawEllipse(int(a.Get(0)),int(a.Get(1)),10,10);
        return;

    }

    // get a second iterator
    typename Container<imfeature>::const_reverse_iterator itb = m_data.rbegin();
    ++itb;
    size_t counter = 0;

    while (itb!=m_data.rend() && counter<length) {

        vec2f a = (*ita).GetLocationAtNativeScale();
        vec2f b = (*itb).GetLocationAtNativeScale();

        p.drawLine(int(a.Get(0)),int(a.Get(1)),int(b.Get(0)),int(b.Get(1)));

        ++ita;
        ++itb;
        ++counter;

    }

}
#endif

template class CTracklet<list>;
template class CTracklet<vector>;
template class CTracklet<CRingBuffer>;


} // end of namespace



