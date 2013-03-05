/*
 * descriptor.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: jbalzer
 */

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
template class CNeighborhoodDescriptor<CRectangle<double>,mat>;
template class CNeighborhoodDescriptor<CRectangle<double>,CDenseVector<bool> >;


}
