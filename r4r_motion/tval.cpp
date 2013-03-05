/*
 * tval.cpp
 *
 *  Created on: Apr 4, 2012
 *      Author: jbalzer
 */

#include "tval.h"
#include "basic.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace cv;
using namespace std;

namespace R4R {

CTrackerValidation::CTrackerValidation(const char* file, cv::Size size):
	m_input(file),
	m_size(size),
	m_fps(m_input.get(CV_CAP_PROP_FPS)) {

}

CTrackerValidation::CTrackerValidation(const char* file, double fps):
	m_input(file),
	m_size(m_input.get(CV_CAP_PROP_FRAME_WIDTH),m_input.get(CV_CAP_PROP_FRAME_HEIGHT)),
	m_fps(fps) {

}

CTrackerValidation::CTrackerValidation(const char* file):
	m_input(file),
	m_size(m_input.get(CV_CAP_PROP_FRAME_WIDTH),m_input.get(CV_CAP_PROP_FRAME_HEIGHT)),
	m_fps(m_input.get(CV_CAP_PROP_FPS)) {

}

bool CTrackerValidation::WriteVideo(const char* file, size_t T) {

	  VideoWriter output(file,CV_FOURCC('D','I','V','X'),m_fps,m_size,true);

	  if(!output.isOpened())
		  return 1;

	  Mat img;

	  for(size_t t=0; t<T; t++) {

		  if(!m_input.grab())
				break;

			m_input.retrieve(img);

			output << img;

	  }

	   return 0;

}

bool CTrackerValidation::WriteTracks(const char* file, std::vector<CTracker*> trackers, size_t tmax) {

	  VideoWriter output(file,CV_FOURCC('M','P','E','G'),m_fps,m_size,true);

	  if(!output.isOpened())
		  return 1;

	  // create an iterator for each tracker
	  vector<CTracker::iterator> tit;
	  for(size_t i=0; i<trackers.size(); i++) {

		  CTracker::iterator it = CTracker::iterator(trackers[i]);
		  tit.push_back(it);

	  }

	  Mat img;
	  size_t t=0;

	  while(t<tmax) {

		  if(!m_input.grab())
				break;

		  m_input.retrieve(img);

		  for(size_t i=0; i<tit.size(); i++) {

			  // check whether tracker is still alive
			  if(tit[i]()) {

				  // current features of the tracker
                  map<shared_ptr<CTracklet>,list<CFeature>::iterator > features = *tit[i];

				  // draw features
                  map<shared_ptr<CTracklet>,list<CFeature>::iterator >::iterator fit;
				  for(fit=features.begin();fit!=features.end();fit++) {
                      (*(fit->second)).Draw(img);

					  // draw tails
				      vec a, b;
                      a = (*(fit->second)).GetLocationAtNativeScale();

					  // iterate backwards
					  while(fit->second!=fit->first->begin()) {

						  fit->second--;

                          b = (*(fit->second)).GetLocationAtNativeScale();

                          line(img,Point2f(a(0),a(1)),Point2f(b(0),b(1)),CFeature::COLORS[(*(fit->second)).GetScale()],2);

						  a = b;

					  }

				  }

				  // increment tracker iterators
				  ++tit[i];

			  }

		  }

		  // write image to output stream
		  output << img;

		  t++;

	  }

	   return 0;

}

void CTrackerValidation::ShowStatistics(std::vector<CTracker*> trackers) {

	for (size_t i=0; i<trackers.size(); i++) {

		cout << "Tracker " << i << ":" << endl;

		size_t sumT(0), maxT(0), noTracks(0);

		for (size_t s = 0; s<trackers[i]->size(); s++) {

			noTracks += trackers[i]->at(s).size();

			list<shared_ptr<CTracklet> >::iterator it;

			for(it=trackers[i]->at(s).begin(); it!=trackers[i]->at(s).end(); it++) {

				if((*it)->size()>maxT)
					maxT = (*it)->size();

				sumT += (*it)->size();

			}

		}

		cout << "Number of tracks: " << noTracks << endl;
		cout << "Maximum life time: " << maxT << " frames" << endl;
		cout << "Average life time: " << ((float)sumT)/((float)noTracks) << " frames" << endl;


	}

}

bool CTrackerValidation::WriteTrack(const char* dir, const char* prefix, shared_ptr<CTracklet> tracklet, size_t hsize) {

	Mat img;

    list<CFeature>::iterator it  = tracklet->begin();

	size_t t = 0;

	while(true) {

		// check whether we can get the image
		if(!m_input.grab() || it==tracklet->end())
			break;

		m_input.retrieve(img);

		if(t>=tracklet->GetCreationTime()) {

			stringstream time;
			time.fill('0');
			time.width(3);
			time << t;


			stringstream filename;

			filename << dir;
			filename << prefix << "_";
			filename << tracklet->GetHash() << "_";
			filename << time.str().c_str() << ".png";

            vec p = (*it).GetLocation();
            double p2s = pow(2,(*it).GetScale());

			Rect roi = Rect(p2s*p(0)-hsize,p2s*p(1)-hsize,2*hsize+1,2*hsize+1);

			Mat subimg = img(roi);

			imwrite(filename.str().c_str(),subimg);

			it++;

		}

		t++;

	}

	return 0;

}


bool CTrackerValidation::SaveToFile(const char* filename, map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > > graph) {

	cout << "Saving covisibility graph to file..." << endl;

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file.\n";
		return 1;

	 }

	map<shared_ptr<CTracklet>,list<shared_ptr<CTracklet> > >::iterator it;
	list<shared_ptr<CTracklet> >::iterator it2;

	for(it=graph.begin(); it!=graph.end(); it++) {

		out << it->first->GetHash() << " ";

		for(it2=it->second.begin(); it2!=it->second.end(); it2++)
			out << (*it2)->GetHash() << " ";

		out << endl;

	}

	out.close();

	return 0;

}

}
