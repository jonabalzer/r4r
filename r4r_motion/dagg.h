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

#ifndef R4RDAGG_H
#define R4RDAGG_H

#include "splinecurve.h"
#include "tracklet.h"

namespace R4R {

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CDescriptorAggregator {

public:

    //! Triggers aggregation.
    virtual void Aggregate(const CTracklet<Container>& tracklet, const string& name, list<imfeature>& result) const = 0;

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CInitFrameAggregator: public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CInitFrameAggregator() {}

    //! \copydoc CDescriptorAggregator<Array,Container>::Aggregate(const CTracklet<Container>&,const string&,list<imfeature>&) const
    void Aggregate(const CTracklet<Container>& tracklet, const string& name, list<imfeature>& result) const;

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CMeanAggregator: public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CMeanAggregator() {}

    //! \copydoc CDescriptorAggregator<Array,Container>::Aggregate(const CTracklet<Container>&,const string&,list<imfeature>&) const
    void Aggregate(const CTracklet<Container>& tracklet, const string& name, list<imfeature>& result) const;

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CSplineInterpolationAggregator: public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CSplineInterpolationAggregator(size_t n, u_int p):m_n(n),m_p(p) {}

    //! \copydoc CDescriptorAggregator<Array,Container>::Aggregate(const CTracklet<Container>&,const string&,list<imfeature>&) const
    void Aggregate(const CTracklet<Container>& tracklet, const string& name, list<imfeature>& result) const;

private:

    size_t m_n;                                         //!< number of control points
    u_int m_p;                                          //!< polynomial degree

};

}
#endif // DAGG_H
