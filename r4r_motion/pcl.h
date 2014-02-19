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

#ifndef R4RPCL_H
#define R4RPCL_H

#include <stdlib.h>

#include "feature.h"
#include "bbox.h"

namespace R4R {

/*! \brief interface for a point cloud or map
 *
 * The class features three template parameters: 1. The container which holds the interest points. In the
 * simplest case, this could be a STL container such as a list. In a more advanced scenario, this could
 * be a tree data structure to facilitate fast lookup, or even something like an (a,b)-tree to store and
 * search very large maps efficiently. 2. The precision for of location coordinates. 3. The dimension of
 * the space. This parameter allows using this class not only as a 3-d map but for example as a collection
 * of interest points in the image plane.
 *
 * The interface should provide the following functionality:
 * - point insertion
 * - closest point search, returns iterator (in terms of geometry)
 * - closest point search (in terms of attached descriptors)
 * - hybrid version of above
 * - maybe point deletion, careful with pointers which are returned upon insertion
 * - read-access to the container
 * - file I/O
 * - calculate bounding box
 *
 */
template<template<class U, class Allocator = std::allocator<U> > class Container,typename T,u_int n=3>
class CPointCloud {

public:

    //! Constructor.
    CPointCloud():m_data() {}

    /*! Insert a new point.
     *
     * \returns pointer to the newly created element in the map
     *
     */
    typename Container<CInterestPoint<T,n> >::const_iterator Insert(const CInterestPoint<T,n>& x) { m_data.push_back(x); return --m_data.end(); }

    //! Read-only access to the data.
    const Container<CInterestPoint<T,n> >& GetData() const { return m_data; }

    //! Write to disk.
    bool WriteToFile(const char* filename) const { return CInterestPoint<T,n>::template SaveToFile<Container>(filename,m_data); }

    //! Finds the geometrically closest point.
    typename Container<CInterestPoint<T,n> >::iterator FindClosestPoint(const CInterestPoint<T,n>& x);

    //! Finds the geometrically closest point.
    typename Container<CInterestPoint<T,n> >::iterator FindClosestPoint(const CInterestPoint<T,n>& x, const std::string& desc);

    //! Computes the boundig box.
    CBoundingBox<T> BoundingBox() const;

    //! Number of points in the cloud.
    size_t Size() const { return m_data.size(); }

    //! Barycenter.
    CVector<T,n> Barycenter() const;

private:

    Container<CInterestPoint<T,n> > m_data;         //!< container holding the data

};



} // end of namespace

#endif
