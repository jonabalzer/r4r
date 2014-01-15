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
#include <set>
#include <typeinfo>

using namespace std;
using namespace cv;

namespace R4R {

template<template<class T, class Allocator = std::allocator<T> > class Container>
CTracker<Container>::CTracker():
    m_data(),
    m_params(nullptr),
    m_global_t(0),
    m_n_active_tracks(0) {}

template<template<class T, class Allocator = std::allocator<T> > class Container>
CTracker<Container>::CTracker(CParameters* params):
    m_data(),
    m_params(params),
    m_global_t(0),
    m_n_active_tracks(0) {}


template<template<class T, class Allocator = std::allocator<T> > class Container>
CTracker<Container>::~CTracker() {

    typename Container<shared_ptr<CTracklet> >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++)
        it->reset();

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracker<Container>::DeleteInvalidTracks() {

    // make a copy of the valid ones, linear complexity
    Container<shared_ptr<CTracklet> > valid;

    typename Container<shared_ptr<CTracklet> >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus())
            valid.push_back((*it));
        else // reset shared ptr
            it->reset();

    }

    // replace data
    m_data = valid;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
shared_ptr<CTracklet> CTracker<Container>::AddTracklet(const imfeature& x) {

    shared_ptr<CTracklet> tracklet(new CTracklet(m_global_t,x));

    m_data.push_back(tracklet);

	return tracklet;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracker<Container>::AddTracklet(CTracklet* tracklet) {

    m_data.push_back(shared_ptr<CTracklet>(tracklet));

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
size_t CTracker<Container>::ActiveCapacity() const {

	size_t result = 0;

    typename Container<shared_ptr<CTracklet> >::const_iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus())
            result++;

    }

	return result;

}

#ifdef QT_GUI_LIB

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracker<Container>::Draw(QImage& img, size_t length) const {

    typename Container<shared_ptr<CTracklet> >::const_iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetStatus()) {

            (*it)->Draw(img, length);

        }

    }

}

#endif

template<template<class T, class Allocator = std::allocator<T> > class Container>
bool CTracker<Container>::SaveToFile(const char* dir, const char* prefix) {

    typename Container<shared_ptr<CTracklet> >::iterator it;

	size_t counter = 0;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        stringstream no;
        no.fill('0');
        no.width(4);
        no << counter;
        counter++;

        stringstream filename;
        filename << string(dir) << string(prefix) << no.str() << ".dat";

        // compress tracklet and reverse order
        (*it)->CompressAndReverse();

        // get access to data
        const vector<imfeature>& data = (*it)->GetData();

        // save to file
        imfeature::SaveToFile<vector>(filename.str().c_str(),data);

    }

	return 0;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
shared_ptr<CTracklet> CTracker<Container>::SearchTracklet(string hash) {

    typename Container<shared_ptr<CTracklet> >::iterator it;

	shared_ptr<CTracklet> result;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        if((*it)->GetHash()==hash)
            result = *it;

    }

	return result;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
shared_ptr<CTracklet> CTracker<Container>::SearchFittestTracklet() {

    typename Container<shared_ptr<CTracklet> >::iterator it;

	size_t max = 0;

	shared_ptr<CTracklet> result;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        size_t lifetime = m_global_t - (*it)->GetCreationTime();

        if(lifetime>max) {

            result = *it;
            max = lifetime;

        }

    }

	return result;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > CTracker<Container>::ComputeCovisibilityGraph() {


	map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > graph;

    typename Container<shared_ptr<CTracklet> >::iterator it, it2;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        list<shared_ptr<CTracklet> > adjlist;
        adjlist.push_back(*it);

        // insert an empty list of neighbors for each tracklet
        graph.insert(pair<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > >(*it,adjlist));

    }

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        for(it2=m_data.begin(); it2!=m_data.end(); it2++) {

            size_t t01 = (*it)->GetCreationTime();
            size_t t02 = (*it2)->GetCreationTime();

            if((*it)!=(*it2) && (t01+(*it)->Size()-1>=t02 || t01<=t02+(*it2)->Size()-1))
                graph[*it].push_back(*it2);

        }

    }

	return graph;

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
void CTracker<Container>::SetAllTrackletsActive() {

    typename Container<shared_ptr<CTracklet> >::iterator it;

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        (*it)->SetStatus(true);

    }

}

template<template<class T, class Allocator = std::allocator<T> > class Container>
vector<size_t> CTracker<Container>::ComputeFeatureDensity(vector<CIntegralImage<size_t> >& imgs) {

    typename Container<shared_ptr<CTracklet> >::iterator it;
    vector<size_t> n(imgs.size());

    for(it=m_data.begin(); it!=m_data.end(); it++) {

        // only consider active tracklets
        if((*it)->GetStatus()) {

            // get the current feature
            imfeature f = (*it)->GetLatestState();

            // see what scale it is at (possibly with rounding)
            u_int s = (u_int)(f.GetScale()+0.5);

            // count
            n[s]++;

            // set Dirac impulse at location in integral image
            vec2f x = f.GetLocation();

            // make sure we have an image at that scale
            if(s<imgs.size())
                imgs[s].AddDensityFast(x.Get(0),x.Get(1),1.0);

        }

    }

    return n;

}

template class CTracker<list>;
template class CTracker<vector>;


} // end of namespace
