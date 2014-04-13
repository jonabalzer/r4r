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

#include "tracker.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;

namespace R4R {

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
CTracker<TrackerContainer,TrackletContainer>::CTracker():
    m_data(),
    m_params(nullptr),
    m_global_t(0) {}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
CTracker<TrackerContainer,TrackletContainer>::CTracker(const CParameters *params):
    m_data(),
    m_params(params),
    m_global_t(0) {}


template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
CTracker<TrackerContainer,TrackletContainer>::~CTracker() {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++)
        it->reset();

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
void CTracker<TrackerContainer,TrackletContainer>::DeleteInvalidTracks() {

    // make a copy of the valid ones, linear complexity
    TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > > valid;

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus())
            valid.push_back((*it));
        else // reset shared ptr
            it->reset();

    }

    // replace data
    m_data = valid;

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
shared_ptr<CTracklet<TrackletContainer> > CTracker<TrackerContainer,TrackletContainer>::AddTracklet(const imfeature& x) {

    shared_ptr<CTracklet<TrackletContainer> > tracklet(new CTracklet<TrackletContainer>(m_global_t,x));

    m_data.push_back(tracklet);

	return tracklet;

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
void CTracker<TrackerContainer,TrackletContainer>::AddTracklet(CTracklet<TrackletContainer>* tracklet) {

    m_data.push_back(shared_ptr<CTracklet<TrackletContainer> >(tracklet));

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
size_t CTracker<TrackerContainer,TrackletContainer>::ActiveCapacity() const {

	size_t result = 0;

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus())
            result++;

    }

	return result;

}

#ifdef QT_GUI_LIB

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
void CTracker<TrackerContainer,TrackletContainer>::Draw(QImage& img, size_t length) const {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus()) {

            if(m_global_t-(*it)->GetCreationTime()>=length)
                (*it)->Draw(img, length);
            else
                (*it)->Draw(img,0);

        }

    }

}

#endif

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
bool CTracker<TrackerContainer,TrackletContainer>::SaveToFile(const char* dir, const char* prefix) {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::iterator it;

	size_t counter = 0;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        stringstream no;
        no.fill('0');
        no.width(4);
        no << counter;
        counter++;

        stringstream filename;
        filename << string(dir) << string(prefix) << no.str() << ".dat";

        // get access to data
        const TrackletContainer<imfeature>& data = (*it)->GetData();

        // save to file
        imfeature::SaveToFile<TrackletContainer>(filename.str().c_str(),data);

    }

	return 0;

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
template<class Array>
list<imfeature> CTracker<TrackerContainer, TrackletContainer>::Aggregate(const CDescriptorAggregator<Array,TrackletContainer>& aggregator, const string& name) const {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it;
    list<imfeature> result;

    for(it=m_data.begin(); it!=m_data.end(); ++it)
        aggregator.Aggregate(**it,name,result);

    return result;

}

template list<imfeature> CTracker<list,CRingBuffer>::Aggregate<matf>(const CDescriptorAggregator<matf,CRingBuffer>& aggregator, const string& name) const;
template list<imfeature> CTracker<list,CRingBuffer>::Aggregate<vecf>(const CDescriptorAggregator<vecf,CRingBuffer>& aggregator, const string& name) const;
template list<imfeature> CTracker<list,list>::Aggregate<matf>(const CDescriptorAggregator<matf,list>& aggregator, const string& name) const;
template list<imfeature> CTracker<list,list>::Aggregate<vecf>(const CDescriptorAggregator<vecf,list>& aggregator, const string& name) const;


template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
shared_ptr<CTracklet<TrackletContainer> > CTracker<TrackerContainer,TrackletContainer>::SearchTracklet(const std::string& hash) const {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it;

    shared_ptr<CTracklet<TrackletContainer> > result;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetHash()==hash)
            result = *it;

    }

	return result;

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
shared_ptr<CTracklet<TrackletContainer> > CTracker<TrackerContainer,TrackletContainer>::SearchFittestTracklet() const {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it;

	size_t max = 0;

    shared_ptr<CTracklet<TrackletContainer> > result;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        size_t lifetime = m_global_t - (*it)->GetCreationTime();

        if(lifetime>=max) {

            result = *it;
            max = lifetime;

        }

    }

	return result;

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
map<shared_ptr<CTracklet<TrackletContainer> >,list<shared_ptr<CTracklet<TrackletContainer> > > > CTracker<TrackerContainer,TrackletContainer>::ComputeCovisibilityGraph() const {

    map<shared_ptr<CTracklet<TrackletContainer> >,list<shared_ptr<CTracklet<TrackletContainer> > > > graph;

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it, it2;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        list<shared_ptr<CTracklet<TrackletContainer> > > adjlist;
        adjlist.push_back(*it);

        // insert an empty list of neighbors for each tracklet
        graph.insert(pair<shared_ptr<CTracklet<TrackletContainer> >,list<shared_ptr<CTracklet<TrackletContainer> > > >(*it,adjlist));

    }

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        for(it2=m_data.begin(); it2!=m_data.end(); it2++) {

            size_t t01 = (*it)->GetCreationTime();
            size_t t02 = (*it2)->GetCreationTime();

            if((*it)!=(*it2) && (t01+(*it)->GetLifetime()-1>=t02 || t01<=t02+(*it2)->GetLifetime()-1))
                graph[*it].push_back(*it2);

        }

    }

	return graph;

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
void CTracker<TrackerContainer,TrackletContainer>::SetAllTrackletsActive() {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        (*it)->SetStatus(true);

    }

}

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
vector<size_t> CTracker<TrackerContainer,TrackletContainer>::ComputeFeatureDensity(vector<CIntImage<size_t> >& imgs) const {

    typename TrackerContainer<shared_ptr<CTracklet<TrackletContainer> > >::const_iterator it;
    vector<size_t> n(imgs.size());

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        // only consider active tracklets
        if((*it)->GetStatus()) {

            // get the current feature
            const imfeature& f = (*it)->GetLatestState();

            // see what scale it is at (possibly with rounding)
            u_int s = u_int(f.GetScale());

            // count
            n[s]++;

            // set Dirac impulse at location in integral image
            const vec2f& x = f.GetLocation();

            // make sure we have an image at that scale
            if(s<imgs.size())
                imgs[s].AddMass(x.Get(1),x.Get(0),1.0);

        }

    }

    return n;

}

template class CTracker<list,list>;
template class CTracker<list,CRingBuffer>;
template class CTracker<list,vector>;
template class CTracker<vector,list>;
template class CTracker<vector,CRingBuffer>;
template class CTracker<vector,vector>;


} // end of namespace
