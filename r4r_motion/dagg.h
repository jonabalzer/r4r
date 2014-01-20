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

#include "tracker.h"
#include "splinecurve.h"
#include <list>

namespace R4R {

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
class CDescriptorAggregator {

public:

    //! Constructor.
    CDescriptorAggregator(CTracker<Container>* tracker, const char* name);

    //! Aggregates over the entire tracker.
    void Aggregate();

    //! Access to the aggregates.
    const list<imfeature>& Get() { return m_aggregate; }

protected:

    //! Aggregates a single tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    //! Copies the feature including only the descriptors specified by #m_name.
    imfeature CopyFeature(const imfeature& x) const;

    CTracker<Container>* m_tracker;                //! tracker to aggregate over
    string m_name;                      //! name of the descriptor to aggregate
    list<imfeature> m_aggregate;         //! list of aggregated features

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CInitFrameAggregator: public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CInitFrameAggregator(CTracker<Container>* tracker, const char* name):CDescriptorAggregator<Array,Container>::CDescriptorAggregator(tracker,name){}

private:

    //! Gets the descriptor from the first feature in the tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    using CDescriptorAggregator<Array,Container>::m_name;
    using CDescriptorAggregator<Array,Container>::m_aggregate;

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CSubsampleAggregator:public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CSubsampleAggregator(CTracker<Container>* tracker, const char* name, size_t n);

private:

    //! Subsamples each tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    size_t m_n;                          //! downsampling factor

    using CDescriptorAggregator<Array,Container>::m_name;
    using CDescriptorAggregator<Array,Container>::m_aggregate;

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CMeanAggregator:public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CMeanAggregator(CTracker<Container>* tracker, const char* name):CDescriptorAggregator<Array,Container>::CDescriptorAggregator(tracker,name){}

private:

    //! Subsamples each tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    using CDescriptorAggregator<Array,Container>::m_name;
    using CDescriptorAggregator<Array,Container>::m_aggregate;

};

template<class Array,template<class T, class Allocator = std::allocator<T> > class Container>
class CSplineInterpolationAggregator:public CDescriptorAggregator<Array,Container> {

public:

    //! Constructor.
    CSplineInterpolationAggregator(CTracker<Container>* tracker, const char* name, size_t n, size_t p);

private:

    //! Subsamples each tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    size_t m_n;                                          //!< number of control points
    size_t m_p;                                          //!< polynomial degree

    using CDescriptorAggregator<Array,Container>::m_name;
    using CDescriptorAggregator<Array,Container>::m_aggregate;

};




}
#endif // DAGG_H
