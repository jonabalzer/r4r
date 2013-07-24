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

} // end of namespace
