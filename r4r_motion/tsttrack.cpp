/*
 * tsthash.cpp
 *
 *  Created on: Mar 17, 2012
 *      Author: jbalzer
 */

#include <iostream>
#include <fstream>

#include "tsttrack.h"
#include "basic.h"
#include "feature.h"
#include "lk.h"
#include "interp.h"

using namespace cv;
using namespace std;

namespace R4R {

CTST::CTST(CParameters params):
	CTracker(params),
	m_root(nullptr),
	m_detector(10000000,m_params.GetDoubleParameter("FEATURE_THRESHOLD"),0,2*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")+1) {
//	m_detector(m_params.GetIntParameter("FEATURE_THRESHOLD"),true) {

	// generate sample points for tests performed in BRIEF descriptor
    CBRIEF::GenerateSamplePoints();

}

void CTST::Detect(vector<Mat>& pyramid) {

	// first create tracklet for every feature point fulfilling distance condition
	for(size_t s = 0; s<=(size_t)m_params.GetIntParameter("SCALE"); s++) {

		// compute integral image for feature counting
		CIntegralImage<size_t> cimg = ComputeFeatureDensity(pyramid[s].cols,pyramid[s].rows,s);

		vector<KeyPoint> keypoints;
		m_detector.detect(pyramid[s],keypoints);

		if(keypoints.size()>0) {

			for(size_t i=0; i<keypoints.size(); i++)
				cimg.AddDensityFast(keypoints[i].pt.x,keypoints[i].pt.y,1);

			cimg.Compute();

			double hsize = m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE");

			for(size_t i=0; i<keypoints.size(); i++) {

				if(cimg.EvaluateFast(keypoints[i].pt.x,keypoints[i].pt.y,hsize,hsize)==1) {

					// create feature
                    CFeature x(keypoints[i].pt.x,keypoints[i].pt.y,s,0);

					// compute descriptor
					CRectangle<double> roi(keypoints[i].pt.x,
										   keypoints[i].pt.y,
										   m_params.GetIntParameter("DESCRIPTOR_HSIZE"),
										   m_params.GetIntParameter("DESCRIPTOR_HSIZE"));

                    shared_ptr<CAbstractDescriptor> brief(new CBRIEF(roi));
					brief->Compute(pyramid[s]);
                    x.AttachDescriptor("BRIEF",brief);

#if COMPUTE_ID_TST == 1
					shared_ptr<CDescriptor> id(new CIdentityDescriptor(roi));
					id->Compute(pyramid[s]);
					x->AttachDescriptor("ID",id);
#endif

					// create new tracklet
                    shared_ptr<CSTTracklet> tracklet(new CSTTracklet(m_global_t,x.GetScale(),x));
					AddTracklet(tracklet);

				}

			}

		}


	}

}

bool CTST::Init(vector<Mat>& pyramid) {

	// create virtual feature/tracklet one scale above desired, FIXME: select scale such that there is actually one feature
	size_t s = m_params.GetIntParameter("SCALE") + 1;

	// root node will be center of image
	vec center(2);
	center(0) = 0.5*(pyramid[0].cols-1);
	center(1) = 0.5*(pyramid[0].rows-1);
	center.Scale((1/pow(2,s)));

	// create the root feature
    CFeature x(center(0),center(1),s);

	// create a new tracklet
    shared_ptr<CSTTracklet> tracklet(new CSTTracklet(m_global_t,x.GetScale(),x));
	AddTracklet(tracklet);

	// make virtual node root of the tree
	m_root = tracklet;

	// perform feature detection
	Detect(pyramid);

	// construct initial tree
	ConstructTree();

	return 0;

}


void CTST::ConstructTree() {

	// now build tree
	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size()-2; s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			shared_ptr<CSTTracklet> parent = FindParent((*it));
			parent->m_children.insert(static_pointer_cast<CSTTracklet>((*it)));

		}

	}

	// attribute all node at coarsest actual level to root
	for(it=at(size()-2).begin(); it!=at(size()-2).end(); it++)
		m_root->m_children.insert(static_pointer_cast<CSTTracklet>((*it)));

}



bool CTST::AddTracklets(vector<Mat>& pyramid) {

	// perform feature detection
	//Detect(img);

	// update the tree

	return 0;


}

shared_ptr<CSTTracklet> CTST::FindParent(shared_ptr<CTracklet> tracklet) {

	list<shared_ptr<CTracklet> >::iterator it;

	vec xc = tracklet->GetLatestLocation();
	size_t scale = tracklet->GetScale();

	// only look at coarser scale
	for(size_t s=scale+1; s<size()-1; s++) {

		// adapt child feature to right scale
		xc.Scale(0.5);

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			vec xp = (*it)->GetLatestLocation();

			double dist = (xp - xc).Norm2();

			if(dist<m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"))
				return static_pointer_cast<CSTTracklet>((*it));

		}

	}

	// if nothing was found, parent must be virtual node
	return m_root;

}

void CTST::DrawChildren(cv::Mat& img, shared_ptr<CSTTracklet> node) {

	vec parent = node->GetLatestLocationAtNativeScale();

	set<shared_ptr<CSTTracklet> >::iterator it;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {


		if((*it)->GetStatus()) {

			size_t s = (*it)->GetScale();

            (*it)->GetLatestState().Draw(img,CFeature::COLORS[s],(m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")+1)*pow(2,s));

		}


		//vec child = (*it)->GetLatestLocationAtNativeScale();

		//line(img,Point(parent(0),parent(1)),Point(child(0),child(1)),CFeature::COLORS[(*it)->GetScale()]);

		DrawChildren(img,(*it));

	}

}

bool CTST::SaveTreeToFile(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file " << filename << "." << endl;
		return 1;

	 }

	out << "digraph G {" << endl;

	// save root
	out << "ranksep=\"1.0 equally\"" << endl;
	out << "ROOT" << " [label=\"\",shape=circle,fillcolor=black,height=0.5,width=0.5];" << endl;

	// save children and the edges recursively
	WriteChildren(out,m_root);

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size()-1; s++) {

		out << "subgraph {" << endl;
		out << "rank = same; ";

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus())
				out << (*it)->GetHash() << "; ";

		}

		out << "}" << endl;

	}


	out << "subgraph {" << endl;
	out << "rank = same; ROOT;" << endl;
	out << "}" << endl;


	out << "}";

	out.close();

	return 0;

}


void CTST::WriteChildren(std::ofstream& os, shared_ptr<CSTTracklet> node) {

	string parent;

    if(!(*node!=*m_root))
		parent = "ROOT";
	else
		parent = node->GetHash();

	set<shared_ptr<CSTTracklet> >::iterator it;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		string child = (*it)->GetHash();

		// first nodes
		//os << child << " [label=\"" << child << "\",shape=circle,fillcolor=black,height=0.5,width=0.5];" << endl;
		os << child << " [label=\"\",shape=circle,fillcolor=black,height=0.25,width=0.25];" << endl;


		// then edges
		os << parent << " -> " << child << endl;

		WriteChildren(os,(*it));

	}

}

bool CTST::RemoveSubTree(shared_ptr<CSTTracklet> node) {

	// check if node is a leaf
	if(node->m_children.size()==0)
		return 0;

	// make a copy of the children
	set<shared_ptr<CSTTracklet> > children = node->m_children;

	// clear children
	node->m_children.clear();

	// recurse with copy of children
	set<shared_ptr<CSTTracklet> >::iterator it;

	for(it=children.begin(); it!=children.end(); it++) {

		// make sure all removed node are also set inactive
		(*it)->SetStatus(false);

		RemoveSubTree(*it);

	}

	return 0;

}

void CTST::PruneTree(shared_ptr<CSTTracklet> node) {

	set<shared_ptr<CSTTracklet> >::iterator it;

	// create container for good children
	set<shared_ptr<CSTTracklet> > children;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		if((*it)->GetStatus()) {

			PruneTree(*it);
			children.insert(*it);

		}
		else
			RemoveSubTree(*it);

	}

	// replace the children container to make sure the bad childs have to be removed themselves
	node->m_children = children;

}


void CTST::TrackSubTree(std::vector<cv::Mat>& pyramid0, vector<Mat>& pyramid1, shared_ptr<CSTTracklet> node, vec t0) {

	size_t hdist = m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE");

	vec u0 = node->GetLatestLocation();

	// predict position
	vec u1 = u0 + t0;

	if(node->GetStatus()) {

	// extract search area
	Rect roid(u1(0)-2*hdist,u1(1)-2*hdist,4*hdist+1,4*hdist+1);				// pad region to have enough data for the detector to function

	// check if we do not exceed image dimensions
	CImageInterpolation::ProjectRect(roid.x,roid.y,roid.width,roid.height,pyramid1[node->GetScale()].cols,pyramid1[node->GetScale()].rows);

	// FIXME: roi extraction 1- or 0-based??
	Mat subimg = pyramid1[node->GetScale()](roid);

	// detect features in search area
	vector<KeyPoint> keypoints;
	m_detector.detect(subimg,keypoints);

/*	if(keypoints.size()==0) {
	namedWindow("Roi");
	imshow("Roi",subimg);
	waitKey(0);

	destroyWindow("Roi");

	}*/


	// just consider those detections in the original roi
	vector<KeyPoint> candidates;
	for(size_t i=0; i<keypoints.size(); i++) {

		if(keypoints[i].pt.x>=hdist && keypoints[i].pt.x<=3*hdist && keypoints[i].pt.y>=hdist && keypoints[i].pt.y<=3*hdist)
			candidates.push_back(keypoints[i]);

	}




	cout << "no of keypoints: " << candidates.size() << endl;

	// check if there is only a single detection
	if(candidates.size()==1) {

		u1(0) = roid.tl().x + candidates[0].pt.x;			// ok, keypoints are 0-based
		u1(1) = roid.tl().y + candidates[0].pt.y;

        CFeature x(u1(0),u1(1),node->GetScale(),0);
		node->Update(x);




	}
	else {

		// kill the node
		node->SetStatus(false);

		//return;

	}

	}

	// motion estimation for children
	set<shared_ptr<CSTTracklet> >::iterator it;
	vec t1 = u1 - u0;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		size_t sd = node->GetScale() - (*it)->GetScale();

		// convert to children's scale
		vec t1s = t1*pow(2,sd);

		if(!node->GetStatus())		// only push down initialization if feature is valid
			t1s.Scale(0);

		TrackSubTree(pyramid0,pyramid1,(*it),t1s);

	}



	//cout << "No of keypoints: " << keypoints.size() << endl;

	// if there are no detections, do nothing (enlarge search region)
/*
	if(keypoints.size()==0) {

		node->SetStatus(false);
		return;

	}

	// acess to last descriptor
	shared_ptr<CDescriptor> desc0 = node->back()->GetDescriptor("BRIEF");
    shared_ptr<CBRIEF> briefdesc0 = static_pointer_cast<CBRIEF>(desc0);

	// array of pointers to tentative descriptors
    CBRIEF* descriptors[keypoints.size()];

	double minhdist = 256;
	size_t imax = 0;

	for(size_t i=0; i<keypoints.size(); i++) {

		// compute BRIEF descriptors for putative feature points
		CRectangle<double> droi(roi.tl().x+keypoints[i].pt.x,roi.tl().y+keypoints[i].pt.y,
								m_params.GetIntParameter("DESCRIPTOR_HSIZE"),
								m_params.GetIntParameter("DESCRIPTOR_HSIZE"));


        descriptors[i] = new CBRIEF(droi,7);
		descriptors[i] ->Compute(pyramid1[node->GetScale()]);

		double hamming = briefdesc0->Distance(*descriptors[i] );

		if(hamming<minhdist) {

			minhdist = hamming;
			imax = i;

		}

	}

	for(size_t i=0; i<keypoints.size(); i++) {

		// delete all suboptimal descriptors
		if(i!=imax) {

			delete descriptors[i];
			descriptors[i] = nullptr;

		}
		else {

			// if the maximum is too big, also delete the minimizer

			if(minhdist>64) {

				delete descriptors[imax];
				node->SetStatus(false);

				return;


			}
			else {	// add feature/descriptor

				u1(0) = roi.tl().x + keypoints[0].pt.x;
				u1(1) = roi.tl().y + keypoints[0].pt.y;

				// create new feature
                CFeature x(u1(0),u1(1),node->GetScale(),minhdist);
				shared_ptr<CDescriptor> brief(descriptors[imax]);
				x->AttachDescriptor("BRIEF",brief);


#if COMPUTE_ID_TST == 1

				CRectangle<double> droi(u1(0),u1(1),
										m_params.GetIntParameter("DESCRIPTOR_HSIZE"),
										m_params.GetIntParameter("DESCRIPTOR_HSIZE"));

				shared_ptr<CDescriptor> id(new CIdentityDescriptor(droi));
				id->Compute(pyramid1[node->GetScale()]);
				x->AttachDescriptor("ID",id);
#endif

				node->Update(x);

				// motion estimation for children
				set<shared_ptr<CSTTracklet> >::iterator it;
				vec t1 = u1 - u0;

				for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

					size_t sd = node->GetScale() - (*it)->GetScale();

					// convert to children's scale
					vec t1s = t1*pow(2,sd);

					if((*it)->GetStatus())
						TrackSubTree(pyramid0,pyramid1,(*it),t1s);

				}


			}

		}


	}
*/

}

bool CTST::Update(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

	set<shared_ptr<CSTTracklet> >::iterator it;

	for(it=m_root->m_children.begin(); it!=m_root->m_children.end(); it++) {

		//if((*it)->GetStatus()) {

			vec t(2);
			t.Scale(0);

			TrackSubTree(pyramid0,pyramid1,(*it),t);

		//}

	}

	m_global_t++;

	return 0;

}


CTSTNode::CTSTNode():
	m_tracklet(nullptr),
	m_parent(nullptr),
	m_children(0) {

}

CTSTNode::CTSTNode(shared_ptr<CTracklet> tracklet):
	m_tracklet(tracklet),
	m_parent(nullptr),
	m_children(0) {

}


CTSTNode::~CTSTNode() {

	m_tracklet->SetStatus(false);
	m_tracklet.reset();

	list<CTSTNode*>::iterator it;

	for(it=m_children.begin(); it!=m_children.end(); it++)
		delete *it;

	m_children.clear();

}

CTSTLK::CTSTLK(CParameters params):
	CTracker(params),
	m_root(nullptr),
//	m_detector(10000000,m_params.GetDoubleParameter("FEATURE_THRESHOLD"),0) {
	m_detector(m_params.GetIntParameter("FEATURE_THRESHOLD")) {

}

bool CTSTLK::Init(vector<Mat>& pyramid) {

	// get scale
	size_t s = m_params.GetIntParameter("SCALE");

	// root node will be center of image
	Point2f center(0.5*(pyramid[0].cols-1),0.5*(pyramid[0].rows-1));
	center = (1/pow(2,s))*center;

	// create the root feature
    CFeature x0(center.x,center.y,s);

	// create a new tracklet
	shared_ptr<CTracklet> tracklet = AddTracklet(x0);

	// create root node
	m_root = new CTSTNode(tracklet);

	// if there is more than one level, construct tree recursively
	if(s>0)
		return Reproduce(pyramid,m_root);

	return 0;

}

bool CTSTLK::ValidateDescendants(Mat& img, CTSTNode* node) {

	// decrease scale
	size_t s = node->m_tracklet->GetScale() - 1;

	// get the current feature state from tracklet attached to the node
	vec xp = node->m_tracklet->GetLatestLocation();

	// get appropriately scaled bounding box around feature
	Rect_<float> roi;

	if(s==(size_t)m_params.GetIntParameter("SCALE")-1) {

		roi = Rect_<float>(m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   img.size().width/pow(2,s)-1-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   img.size().height/pow(2,s)-1-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

	}
	else {

		roi = Rect_<float>(2*(xp(0)-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")),
						   2*(xp(1)-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")),
						   4*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   4*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

	}

	// compute integral image to count the number of features in this neighborhood at next finer resolution
	CIntegralImage<size_t> cimg = CIntegralImage<size_t>(roi.width,roi.height);
	cimg.Clear();

	list<CTSTNode*>::iterator it;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		// get the latest state of the children
		vec xc = (*it)->m_tracklet->GetLatestLocation();

		// set input of integral image to one
		cimg.AddDensityFast(xc(0)-roi.tl().x,xc(1)-roi.tl().y,1);

	}

	// compute the counting image
	cimg.Compute();

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		vec xc = (*it)->m_tracklet->GetLatestLocation();

		size_t no = cimg.EvaluateFast(xc(0)-roi.tl().x,
						 	 	 	  xc(1)-roi.tl().y,
						 	 	 	  m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),										// FIXME: add factor 2 that enforces that the windows are disjoint
						 	 	 	  m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

		if(no!=1)
			(*it)->m_tracklet->SetStatus(false);		// if we find a bad track, mark it
		else
			ValidateDescendants(img,(*it));				// if not, check the next level of nodes

	}

	return 0;

}

bool CTSTLK::Update(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

	// collects all nodes in the tree
	vector<list<CTSTNode*> > nodes(m_params.GetIntParameter("SCALE")+1);
	CollectSiblings(nodes,m_root);

	if(nodes.size()==0)
		return 0;

	// start just below the root node
	for(int s=m_params.GetIntParameter("SCALE")-1; s>=0; s--) {

		// check whether there is something to track at scale s
		if(nodes[s].size()>0) {

			// preallocate point arrays to make this faster
			vector<Point2f> points0(nodes[s].size());
			vector<Point2f> points1(nodes[s].size());
			vector<uchar> status(nodes[s].size());
			vector<float> error(nodes[s].size());

			list<CTSTNode*>::iterator it;
			size_t counter = 0;

			for(it=nodes[s].begin(); it!=nodes[s].end(); it++) {

				vec x0 = (*it)->m_tracklet->GetLatestLocation();
				points0[counter] = Point2f(x0(0),x0(1));

				// do not use initial value of invalid parent
				if(!(*it)->m_parent->m_tracklet->GetStatus())
					points1[counter] = points0.back();
				else {

					// get velocity from parent which has been already updated
					vec vp = (*it)->m_parent->m_tracklet->GetLatestVelocity();

					// get current location of node at current level
					vec x = (*it)->m_tracklet->GetLatestLocation();

					// predict the new state of node at current level by parent's velocity scaled correctly
					points1[counter] = Point2f(x(0)+2*vp(0),x(1)+2*vp(0));

				}

				counter++;

			}

			// perform the actual tracking
			calcOpticalFlowPyrLK(pyramid0[s],
								 pyramid1[s],
								 points0,
								 points1,
								 status,
								 error,
								 Size(2*m_params.GetIntParameter("TRACKING_HSIZE")+1,2*m_params.GetIntParameter("TRACKING_HSIZE")+1),
								 m_params.GetIntParameter("LK_PYRAMID_LEVEL"),
								 TermCriteria(TermCriteria::COUNT|TermCriteria::EPS,
										 	  m_params.GetIntParameter("MAX_ITER"),
										 	  m_params.GetDoubleParameter("ACCURACY")),
							     m_params.GetDoubleParameter("LAMBDA"),
								 OPTFLOW_USE_INITIAL_FLOW);

			counter = 0;

			// update the tracks
			for(it=nodes[s].begin(); it!=nodes[s].end(); it++) {

				// update only if the low-level tracking result is ok (do not use SSD as a quality measure!)
				if(status[counter] && points1[counter].x>=1 && points1[counter].x<=pyramid0[s].cols-2 && points1[counter].y>=1 && points1[counter].y<=pyramid0[s].rows-2) {

					// create new feature
                    CFeature x(points1[counter].x,points1[counter].y,s,error[counter]);

					// append feature at the end of the track
					(*it)->m_tracklet->Update(x);

				} else
					(*it)->m_tracklet->SetStatus(false);

				counter++;

			}

		}

	}

	m_global_t++;

	return 0;

}

void CTSTLK::CollectSiblings(vector<list<CTSTNode*> >& nodes, CTSTNode* node) {

	nodes[node->m_tracklet->GetScale()].push_back(node);

	list<CTSTNode*>::iterator it;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++)
		CollectSiblings(nodes,(*it));

}

bool CTSTLK::Prune(CTSTNode* node) {

	list<CTSTNode*>::iterator it;
	list<CTSTNode*> children;

	// breadth first delete
	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		if(!(*it)->m_tracklet->GetStatus())
			delete (*it);
		else
			children.push_back((*it));

	}

	// update the children list
	node->m_children = children;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++)
		Prune((*it));

	return 0;

}

void CTSTLK::Clean(vector<Mat>& pyramid0, vector<Mat>& pyramid1) {

	// mark childless children of root for deletion
	list<CTSTNode*>::iterator it;
	for(it=m_root->m_children.begin(); it!=m_root->m_children.end(); it++) {

		if((*it)->m_children.size()==0)
			(*it)->m_tracklet->SetStatus(false);

	}

	// if there are more than one level, recursively check the nodes in the tree
	if(m_params.GetIntParameter("SCALE")>0) {

		// find bad nodes and mark them
		ValidateDescendants(pyramid1[0],m_root);

		// clean up the tree
		Prune(m_root);

	}

	// FIXME: counter number of trees recursively

}

bool CTSTLK::AddTracklets(vector<Mat>& pyramid) {

	// if there are more than one level, recursively expand the tree
	if(m_params.GetIntParameter("SCALE")>0)
		return Reproduce(pyramid,m_root);

	return 0;

}

bool CTSTLK::Reproduce(vector<Mat>& pyramid, CTSTNode* node) {

	// decrease scale
	size_t s = node->m_tracklet->GetScale() - 1;

	// get the current feature state from tracklet attached to the node
	vec xp = node->m_tracklet->GetLatestLocation();

	// get appropriately scaled bounding box around feature
	Rect_<float> roi;

	if(s==(size_t)m_params.GetIntParameter("SCALE")-1) {

		roi = Rect_<float>(m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   pyramid[s].size().width-1-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   pyramid[s].size().height-1-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

	}
	else {

		roi = Rect_<float>(2*(xp(0)-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")),
						   2*(xp(1)-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")),
						   4*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
						   4*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

	}

	// if roi is inside the image, extract subimage at new scale
	Mat subimg;

	if(0<=roi.x && roi.x + roi.width<=pyramid[s].cols && 0<=roi.y && roi.y+roi.height<=pyramid[s].rows)
		subimg = pyramid[s](roi);
	else
		return 1;

	// look for features
	vector<KeyPoint> keypoints;
	m_detector.detect(subimg,keypoints);

	// if there are no more features, node is a leave
	if(keypoints.size()==0)
		return 1;

	// allocate space for integral image
	CIntegralImage<size_t> cimg = CIntegralImage<size_t>(subimg.size().width,subimg.size().height);
	cimg.Clear();

	for(size_t i=0; i<keypoints.size(); i++)
		cimg.AddDensityFast(keypoints[i].pt.x,keypoints[i].pt.y,1);

	// locations of existing points at that scale
	list<CTSTNode*>::iterator it;

	for(it=node->m_children.begin(); it!=node->m_children.end(); it++) {

		// get the latest state of the children w.r.t. bounding box origin
		vec xc = (*it)->m_tracklet->GetLatestLocation();

		cimg.AddDensityFast(xc(0),xc(1),1);

	}

	// compute the integral image
	cimg.Compute();

	for(size_t i=0; i<keypoints.size(); i++) {

		// evaluate integral image, using down-scaled window, newly-detected keypoints are integral
		size_t no = cimg.EvaluateFast(keypoints[i].pt.x,
						 	 	 	  keypoints[i].pt.y,
						 	 	 	  m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),			// FIXME: add factor 2 that enforces that the windows are disjoint
						 	 	 	  m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

		if(no==1) {

			// transform feature location back to image coordinates
			Point2f feature = roi.tl() + keypoints[i].pt;

			// create new feature and a track for it
            CFeature x(feature.x,feature.y,s);
			shared_ptr<CTracklet> tracklet = AddTracklet(x);

			// insert new node into the tree
			CTSTNode* child = new CTSTNode(tracklet);
			child->m_parent = node;
			node->m_children.push_back(child);

			// if there are levels left, continue
			if(s>0)
				Reproduce(pyramid,child);


		}

	}

	return 0;

}

/*

void CTSTLK::Draw(Mat& img) {

	list<shared_ptr<CTracklet> >::iterator it;

	for(size_t s=0; s<size(); s++) {

		for(it=at(s).begin(); it!=at(s).end(); it++) {

			if((*it)->GetStatus() && s!=(size_t)m_params.GetIntParameter("SCALE")) {

				vec x = (*it)->GetLatestLocation();

				Rect_<float> roi(pow(2,s)*(x(0)-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")),
								 pow(2,s)*(x(1)-m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE")),
								 pow(2,s+1)*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"),
								 pow(2,s+1)*m_params.GetIntParameter("MINIMAL_FEATURE_HDISTANCE"));

				rectangle(img,roi.tl(),roi.br(),CFeature::COLORS[s],1);

			}

		}

	}

}
*/





}	// end of namespace

