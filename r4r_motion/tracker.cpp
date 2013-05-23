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
	vector<list<shared_ptr<CTracklet> > >(),
    m_global_t(0),
    m_n_active_tracks(0) {

}

CTracker::CTracker(CParameters* params):
    vector<list<shared_ptr<CTracklet> > >(params->GetIntParameter("SCALE")+1),
    //vector<list<shared_ptr<CTracklet> > >(6), // hard-coded limit
    m_params(params),
	m_global_t(0)
{}

CTracker::~CTracker() {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++)
			it->reset();

	}

	clear();

}

void CTracker::DeleteInvalidTracks() {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		vector<shared_ptr<CTracklet> > todel;

		// mark for deletion
		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if(!(*it)->GetStatus() && (*it)!=nullptr) {

                todel.push_back(move(*it));

			}

		}

		// delete and remove from list
		for(size_t i=0; i<todel.size(); i++) {

			todel[i].reset();
			at(s).remove(todel[i]);

		}

	}

}

shared_ptr<CTracklet> CTracker::AddTracklet(CFeature x) {

    shared_ptr<CTracklet> tracklet(new CTracklet(m_global_t,x.GetScale(),x));

	// discrepancy between size of vector and s
	int d = tracklet->GetScale() - size();

	// insert empty sets
	for(int i=0; i<=d; i++)
		push_back(list<shared_ptr<CTracklet> >());

	// insert tracklet at s
	at(tracklet->GetScale()).push_back(tracklet);

	return tracklet;

}

void CTracker::AddTracklet(std::shared_ptr<CTracklet> tracklet) {

	// discrepancy between size of vector and s
	int d = tracklet->GetScale() - size();

	// insert empty sets
	for(int i=0; i<=d; i++)
		push_back(list<shared_ptr<CTracklet> >());

	// insert tracklet at s
	at(tracklet->GetScale()).push_back(tracklet);

}

void CTracker::Clean(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1) {

    list<shared_ptr<CTracklet> >::iterator it;
    size_t counter = 0;

    for(size_t s=0; s<size(); s++) {

            // mark for deletion
        for(it=at(s).begin(); it!=at(s).end(); it++) {

            if((*it)->GetStatus())
                counter++;

        }

    }

    m_n_active_tracks = counter;

}

size_t CTracker::Capacity() {

	size_t result = 0;

	for(size_t s=0; s<size(); s++)
		result += at(s).size();

	return result;

}

size_t CTracker::ActiveCapacity() {

	size_t result = 0;

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus())
				result++;


		}

	}

	return result;

}

void CTracker::Draw(Mat& img) {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

                CFeature x = (*it)->GetLatestState();

				if(m_global_t - (*it)->GetCreationTime()<10) {

					// draw additional inner circle
                    x.Draw(img,CFeature::COLORS[(*it)->GetScale()],2);
                    x.Draw(img);

				}
				else
                    x.Draw(img);


			}

		}

	}

}

void CTracker::DrawTails(cv::Mat& img, size_t length) {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

				(*it)->Draw(img, length);

			}

		}

	}

}

void CTracker::DeleteDescriptors() {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			(*it)->DeleteDescriptors();

		}

	}

}

bool CTracker::SaveToFile(const char* dir, const char* prefix) {

	cout << "Writing tracks to file..." << endl;

	list<shared_ptr<CTracklet> >::iterator it;

	size_t counter = 0;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			stringstream no;
			no.fill('0');
			no.width(4);
			no << counter;
			counter++;

			stringstream filename;
			filename << dir << prefix << no.str() << ".dat";

            CFeature::SaveToFile(filename.str().c_str(),*(*it));

		}

	}

	return 0;

}

shared_ptr<CTracklet> CTracker::SearchTracklet(CFeature x0, size_t t0) {

	list<shared_ptr<CTracklet> >::iterator it;

	shared_ptr<CTracklet> result;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

            if((*it)->GetCreationTime()==t0 && x0.operator ==((*it)->front()))
				result = *it;

		}

	}

	return result;

}

shared_ptr<CTracklet> CTracker::SearchFittestTracklet() {

	list<shared_ptr<CTracklet> >::iterator it;

	size_t max = 0;

	shared_ptr<CTracklet> result;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->size()>max) {

				result = *it;

				max = result->size();

			}

		}

	}

	return result;

}


map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > CTracker::ComputeCovisibilityGraph() {

	cout << "Computing covisibility graph..." << endl;

	map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > graph;

	list<shared_ptr<CTracklet> >::iterator it, it2;

	// each tracklet is covisible with itself
	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			list<shared_ptr<CTracklet> > adjlist;
			adjlist.push_back(*it);

			graph.insert(pair<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > >(*it,adjlist));

		}

	}

	// now compute general relationships
	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			for(size_t t=0; t<size(); t++) {

				for(it2=at(s).begin(); it2!=at(s).end(); it2++) {

					size_t t01 = (*it)->GetCreationTime();
					size_t t02 = (*it2)->GetCreationTime();

					if((*it)!=(*it2) && (t01+(*it)->size()-1>=t02 || t01<=t02+(*it2)->size()-1))
						graph[*it].push_back(*it2);

				}

			}

		}

	}

	return graph;

}


void CTracker::SetAllTrackletsActive() {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			(*it)->SetStatus(true);

		}

	}


}



CIntegralImage<size_t> CTracker::ComputeFeatureDensity(size_t width, size_t height, size_t s) {

	// init image
	CIntegralImage<size_t> img(width,height);

	// check whether there are tracklets at that scale at all
	if(at(s).size()>0) {

		list<shared_ptr<CTracklet> >::iterator it;

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus()) {

				vec pt = (*it)->GetLatestLocation();

				img.AddDensityFast(pt(0),pt(1),1.0);

			}

		}

	}

	return img;

}


void CTracker::OnMouseSelectBoundingBox(int event, int x, int y, int flags, void* params) {

	int* ploc = (int*)params;

	if (event == CV_EVENT_LBUTTONUP) {

		ploc[0] = x;
		ploc[1] = y;
		ploc[2] = -1;

	}

	if (event == CV_EVENT_RBUTTONUP) {

		ploc[0] = x;
		ploc[1] = y;
		ploc[2] = 1;

	}

}


Rect CTracker::GetManualBoundingBox(Mat& img) {

	namedWindow("Init",CV_WINDOW_AUTOSIZE|CV_GUI_NORMAL);

	// set intial bounding box
	cout << "Select initial bounding box by clicking on top-left/bottom-right corners..." << endl;

	displayOverlay("Init","Select initial bounding box by clicking on top-left/bottom-right corners...");

	int p[3];

	cvSetMouseCallback("Init",CTracker::OnMouseSelectBoundingBox,p);

	Point2i tl(0,0), br(0,0);

	while (1) {

		imshow("Init",img);

		if(p[2]==-1) {

			tl.x = p[0];
			tl.y = p[1];

		}

		if(p[2]==1 && p[0]>tl.x && p[1]>tl.y) {

			br.x = p[0];
			br.y = p[1];
			break;

		}

		waitKey(33);

	}

	rectangle(img,tl,br,Scalar(255,255,255),2);

	imshow("Init",img);

	cout << "Press any key to continue..." << endl;

	waitKey(0);

	destroyWindow("Init");

	return Rect(tl.x,tl.y,br.x-tl.x,br.y-tl.y);

}



CTracker::iterator::iterator(CTracker* tracker):
	m_tracker(tracker),
	m_t(0) {

	AddTracklets();

}


void CTracker::iterator::UpdateTracklets() {

    map<shared_ptr<CTracklet>,list<CFeature>::iterator >::iterator it;

	for(it=m_data.begin(); it!=m_data.end(); it++)
		(it->second)++;

}

void CTracker::iterator::AddTracklets() {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<m_tracker->size(); s++) {

		for(it=m_tracker->at(s).begin(); it!=m_tracker->at(s).end(); it++) {

			if((*it)->GetCreationTime()==m_t) {

                list<CFeature>::iterator itt = it->get()->begin();

                m_data.insert(pair<shared_ptr<CTracklet>,list<CFeature>::iterator >(*it,itt));

			}

		}

	}

}

void CTracker::iterator::Clean() {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<m_tracker->size(); s++) {

		for(it=m_tracker->at(s).begin(); it!=m_tracker->at(s).end(); it++) {

			if(m_t>(*it)->GetCreationTime()+(*it)->size()-1)
				m_data.erase(*it);

		}

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
