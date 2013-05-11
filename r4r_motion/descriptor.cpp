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

#include "descriptor.h"
#include "rect.h"

using namespace cv;
using namespace std;

namespace R4R {


template <class Array>
CDescriptor<Array>::CDescriptor():
    m_container() {}

template <class Array>
CDescriptor<Array>::CDescriptor(const Array& container):
    m_container(container) {}

template class CDescriptor<mat>;
template class CDescriptor<matf>;
template class CDescriptor<vec>;
template class CDescriptor<vecf>;
template class CDescriptor<CDenseArray<int> >;
template class CDescriptor<CDenseArray<bool> >;
template class CDescriptor<CDenseArray<size_t> >;
template class CDescriptor<CDenseVector<bool> >;


template <class Rect, class Array>
CNeighborhoodDescriptor<Rect,Array>::CNeighborhoodDescriptor(const Rect& roi):
        CDescriptor<Array>::CDescriptor(),
		m_roi(roi) {}

template <class Rect, class Array>
CNeighborhoodDescriptor<Rect,Array>::CNeighborhoodDescriptor(const Array& container, const Rect& roi):
    CDescriptor<Array>::CDescriptor(container),
    m_roi(roi) {}

template class CNeighborhoodDescriptor<CRectangle<double>,matf>;
template class CNeighborhoodDescriptor<CRectangle<double>,vecf>;
template class CNeighborhoodDescriptor<CRectangle<double>,mat>;
template class CNeighborhoodDescriptor<CRectangle<double>,CDenseVector<bool> >;


CDescriptorFileHeader::CDescriptorFileHeader():
    m_scale(0),
    m_name(),
    m_type(ETYPE::NA),
    m_comment() {

    m_location[0] = 0;
    m_location[1] = 0;
    m_size[0] = 0;
    m_size[1] = 0;

}

CDescriptorFileHeader::CDescriptorFileHeader(const CFeature& feature):
    m_scale(feature.m_scale),
    m_name(),
    m_type(ETYPE::NA),
    m_comment() {

    m_location[0] = feature.m_location.Get(0);
    m_location[1] = feature.m_location.Get(1);
    m_size[0] = 0;
    m_size[1] = 0;

}

bool CDescriptorFileHeader::SetDescriptor(CFeature& feature, const char* name) {

    // check if descriptor exists
    if(feature.HasDescriptor(name)) {

        shared_ptr<CAbstractDescriptor> pdesc = feature.GetDescriptor(name);
        m_size[0] = pdesc->NRows();
        m_size[1] = pdesc->NCols();

        m_name = string(name);
        m_type = pdesc->GetType();

    }
    else
        return 1;

    return 0;

}

ostream& operator<<(std::ostream& os, CDescriptorFileHeader& x) {

    os << x.m_location[0] << " " << x.m_location[1] << endl;
    os << x.m_scale << endl;
    os << x.m_name << endl;
    os << x.m_size[0] << " " << x.m_size[1] << endl;
    os << (int)x.m_type << endl;
    os << x.m_comment;

    return os;

}

istream& operator>>(std::istream& is, CDescriptorFileHeader& x) {

    is >> x.m_location[0];
    is >> x.m_location[1];
    is.get();

    is >> x.m_scale;
    is.get();

    is >> x.m_name;
    is.get();

    is >> x.m_size[0];
    is >> x.m_size[1];

    int temp;
    is >> temp;
    x.m_type = (ETYPE)temp;
    is.get();

    is >> x.m_comment;

    return is;

}

}
