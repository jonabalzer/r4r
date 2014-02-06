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

#include <list>
#include <limits>

#include "pcl.h"

using namespace std;

namespace R4R {

template<template<class U, class Allocator = std::allocator<U> > class Container,typename T,u_int n>
typename Container<CInterestPoint<T,n> >::iterator CPointCloud<Container,T,n>::FindClosestPoint(const CInterestPoint<T,n>& x) {



}

template<template<class U, class Allocator = std::allocator<U> > class Container,typename T,u_int n>
typename Container<CInterestPoint<T,n> >::iterator CPointCloud<Container,T,n>::FindClosestPoint(const CInterestPoint<T,n>& x, const string& desc) {



}

template<template<class U, class Allocator = std::allocator<U> > class Container,typename T,u_int n>
CBoundingBox<T> CPointCloud<Container,T,n>::BoundingBox() const {



    typename Container<CInterestPoint<T,n> >::const_iterator it;

    CVector<T,n> lower, upper;

    if(this->Size()==0)
        return CBoundingBox<T>(lower,upper);

    for(u_int i=0; i<n; i++) {

        lower(i) = std::numeric_limits<float>::max();
        upper(i) = -std::numeric_limits<float>::max();

    }

    for(it=m_data.begin(); it!=m_data.end(); ++it) {

        const CVector<T,n>& point = it->GetLocation();

        for(u_int i=0; i<n; i++) {

            if(point.Get(i)>=upper.Get(i))
                upper(i) = point.Get(i);

            if(point.Get(i)<=lower.Get(i))
                lower(i) = point.Get(i);

        }

    }

    return CBoundingBox<T>(lower,upper);

}

template class CPointCloud<list,float,3>;

}
