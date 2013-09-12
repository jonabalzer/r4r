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
#include <GL/gl.h>
#include <GL/glut.h>

using namespace std;
using namespace cv;

namespace R4R {

CTracker::CTracker():
    list<shared_ptr<CTracklet> >(),
    m_global_t(0),
    m_n_active_tracks(0) {}

CTracker::CTracker(CParameters* params):
    list<shared_ptr<CTracklet> >(),
    m_params(params),
    m_global_t(0) {}

CTracker::~CTracker() {

	list<shared_ptr<CTracklet> >::iterator it;

    for(it=begin(); it!=end(); it++)
        it->reset();

	clear();

}

void CTracker::DeleteInvalidTracks() {

    list<shared_ptr<CTracklet> >::iterator it;
    vector<shared_ptr<CTracklet> > todel;

    // mark for deletion
    for(it=begin(); it!=end(); it++) {

         if(!(*it)->GetStatus() && (*it)!=nullptr)
             todel.push_back(move(*it));

    }

    // delete and remove from list
    for(size_t i=0; i<todel.size(); i++) {

        // the tracklet will be destroyed if this is the last reference
        todel[i].reset();
        remove(todel[i]);

    }

}

shared_ptr<CTracklet> CTracker::AddTracklet(imfeature x) {

    shared_ptr<CTracklet> tracklet(new CTracklet(m_global_t,x));

    push_back(tracklet);

	return tracklet;

}

size_t CTracker::ActiveCapacity() const {

	size_t result = 0;

    list<shared_ptr<CTracklet> >::const_iterator it;

    for(it=begin(); it!=end(); it++) {

        if((*it)->GetStatus())
            result++;

    }

	return result;

}

#ifdef QT_GUI_LIB

void CTracker::Draw(QImage& img, size_t length) const {

    list<shared_ptr<CTracklet> >::const_iterator it;

    for(it=begin(); it!=end(); it++) {

        if((*it)->GetStatus()) {

            (*it)->Draw(img, length);

        }

    }

}

#endif

bool CTracker::SaveToFile(const char* dir, const char* prefix) {

	cout << "Writing tracks to file..." << endl;

	list<shared_ptr<CTracklet> >::iterator it;

	size_t counter = 0;

    for(it=begin(); it!=end(); it++) {

        stringstream no;
        no.fill('0');
        no.width(4);
        no << counter;
        counter++;

        stringstream filename;
        filename << string(dir) << string(prefix) << no.str() << ".dat";

        imfeature::SaveToFile(filename.str().c_str(),*(*it));

    }

	return 0;

}

shared_ptr<CTracklet> CTracker::SearchTracklet(imfeature x0, size_t t0) {

	list<shared_ptr<CTracklet> >::iterator it;

	shared_ptr<CTracklet> result;

    for(it=begin(); it!=end(); it++) {

        if((*it)->GetCreationTime()==t0 && (*it)->front()==x0)
            result = *it;

    }

	return result;

}

shared_ptr<CTracklet> CTracker::SearchFittestTracklet() {

	list<shared_ptr<CTracklet> >::iterator it;

	size_t max = 0;

	shared_ptr<CTracklet> result;

    for(it=begin(); it!=end(); it++) {

        if((*it)->size()>max) {

            result = *it;

            max = result->size();

        }

    }

	return result;

}


map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > CTracker::ComputeCovisibilityGraph() {

	cout << "Computing covisibility graph..." << endl;

	map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > graph;

	list<shared_ptr<CTracklet> >::iterator it, it2;

    for(it=begin(); it!=end(); it++) {

        list<shared_ptr<CTracklet> > adjlist;
        adjlist.push_back(*it);

        graph.insert(pair<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > >(*it,adjlist));

    }

    for(it=begin(); it!=end(); it++) {

        for(size_t t=0; t<size(); t++) {

            for(it2=begin(); it2!=end(); it2++) {

                size_t t01 = (*it)->GetCreationTime();
                size_t t02 = (*it2)->GetCreationTime();

                if((*it)!=(*it2) && (t01+(*it)->size()-1>=t02 || t01<=t02+(*it2)->size()-1))
                    graph[*it].push_back(*it2);

            }

        }

    }

	return graph;

}


void CTracker::SetAllTrackletsActive() {

	list<shared_ptr<CTracklet> >::iterator it;

    for(it=begin(); it!=end(); it++) {

        (*it)->SetStatus(true);

    }

}

vector<size_t> CTracker::ComputeFeatureDensity(vector<CIntegralImage<size_t> >& imgs) {

    list<shared_ptr<CTracklet> >::iterator it;
    vector<size_t> n(imgs.size());

    for(it=begin(); it!=end(); it++) {

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

CTracker::iterator::iterator(CTracker* tracker):
	m_tracker(tracker),
	m_t(0) {

	AddTracklets();

}


void CTracker::iterator::UpdateTracklets() {

    map<shared_ptr<CTracklet>,list<imfeature>::iterator >::iterator it;

	for(it=m_data.begin(); it!=m_data.end(); it++)
		(it->second)++;

}

void CTracker::iterator::AddTracklets() {

	list<shared_ptr<CTracklet> >::iterator it;

    for(it=m_tracker->begin(); it!=m_tracker->end(); it++) {

        if((*it)->GetCreationTime()==m_t) {

            list<imfeature>::iterator itt = it->get()->begin();

            m_data.insert(pair<shared_ptr<CTracklet>,list<imfeature>::iterator >(*it,itt));

        }

    }

}

void CTracker::iterator::Clean() {

	list<shared_ptr<CTracklet> >::iterator it;

    for(it=m_tracker->begin(); it!=m_tracker->end(); it++) {

        if(m_t>(*it)->GetCreationTime()+(*it)->size()-1)
            m_data.erase(*it);

    }

}

void CTracker::iterator::operator++() {

	m_t++;

	UpdateTracklets();

	Clean();

	AddTracklets();

}


void CTracker::iterator::Advance(size_t step) {

	size_t k = 0;

	while(k<step) {

		if(this->operator ()())
			this->operator ++();
		else
			break;

		k++;

	}

}




} // end of namespace
