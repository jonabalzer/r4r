/*
 * feature.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: jbalzer
 */

#include "feature.h"
#include "descriptor.h"
#include <stdio.h>
#include <assert.h>

using namespace cv;
using namespace std;

namespace R4R {

const Scalar CFeature::COLORS[12] = { Scalar(0,255,0), Scalar(255,0,0), Scalar(0,0,255),
									 Scalar(255,255,255), Scalar(255,255,0), Scalar(255,0,255),
									 Scalar(0,128,0), Scalar(128,0,0), Scalar(0,0,128),
									 Scalar(128,128,128), Scalar(128,128,0), Scalar(128,0,128) };

CFeature::CFeature():
        m_location(2),
		m_scale(0),
		m_quality(0) {}

CFeature::CFeature(vec loc, size_t scale, double quality):
		m_location(loc),
		m_scale(scale),
		m_quality(quality) {}

CFeature::CFeature(double x, double y, size_t scale, double quality):
		m_location(2),
		m_scale(scale),
		m_quality(quality) {

	m_location(0) = x;
	m_location(1) = y;

}

ostream& operator<<(ostream& os, CFeature& x) {

    os << "Location: ";

    // print location depending on number of dimensions
    for(size_t i=0; i<x.m_location.NElems()-1; i++)
        os << x.m_location.Get(i) << " ";

    os << x.m_location.Get(x.m_location.NElems()-1) << endl;

    // scale
    os << "Scale: " << x.m_scale << endl;

    // quality
    os << "Quality: " << x.m_quality << endl;

    // number of descriptors that follow
    os << "Number of descriptors: " << x.NoDescriptors() << endl;

    return os;

}

ofstream& operator<<(ofstream& os, CFeature& x) {

	// print location depending on number of dimensions
	for(size_t i=0; i<x.m_location.NElems()-1; i++)
		os << x.m_location.Get(i) << " ";

	os << x.m_location.Get(x.m_location.NElems()-1) << endl;

	// scale
	os << x.m_scale << endl;

	// quality
	os << x.m_quality << endl;

	// number of descriptors that follow
	os << x.NoDescriptors() << endl;

    map<string,shared_ptr<CAbstractDescriptor> >::iterator it;

	// print all descriptors
    for(it = x.m_descriptors.begin(); it!=x.m_descriptors.end(); it++) {

		os << it->first.c_str() << endl;
		it->second->Write(os);
		os << endl;

	}

	return os;

}

ifstream& operator>>(std::ifstream& is, CFeature& x) {

	// read location
	is >> x.m_location(0);
	is >> x.m_location(1);
	is.get();

	// scale
	is >> x.m_scale;
	is.get();

	// read quality
	is >> x.m_quality;
	is.get();

	// read number of descriptors
	size_t nod;
	is >> nod;
	is.get();

	for(size_t k=0; k<nod; k++) {

        // name
        string name;
        getline(is,name);

        // read all as floats
        CDescriptor<vecf>* pdesc = new CDescriptor<vecf>();
        pdesc->Read(is);
        x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

	}

	return is;

}


bool CFeature::operator<(CFeature& x) {

	vec pt0 = this->GetLocation();
	vec pt1 = x.GetLocation();

	if(pt1(0)<pt0(0))
		return true;
	else
		return pt1(0)==pt0(0) && pt1(1)<pt0(1);

}


bool CFeature::operator!=(CFeature& x) {

	vec pt0 = this->GetLocation();
	vec pt1 = x.GetLocation();

	return pt0(0)!=pt1(0) && pt0(1)!=pt1(1);

}

CFeature::~CFeature() {

	DeleteDescriptors();

}

void CFeature::DeleteDescriptors() {

	map<string,shared_ptr<CAbstractDescriptor> >::iterator it;

	for(it=m_descriptors.begin(); it!=m_descriptors.end(); it++)
		it->second.reset();

	m_descriptors.clear();


}

void CFeature::Draw(Mat& img, Scalar color, size_t size) {

	if(m_location.NElems()==2)
		circle(img,Point(pow(2,(double)m_scale)*m_location.Get(0),pow(2,(double)m_scale)*m_location.Get(1)),size,color);

}

void CFeature::Draw(cv::Mat& img) {

    if(m_location.NElems()==2)
		Draw(img,COLORS[m_scale],3*pow(2,m_scale));

}

void CFeature::AttachDescriptor(const char* id, shared_ptr<CAbstractDescriptor> descriptor) {

	m_descriptors.insert(pair<string,shared_ptr<CAbstractDescriptor> >(id,descriptor));

}

shared_ptr<CAbstractDescriptor> CFeature::GetDescriptor(const char* name) {

	if(m_descriptors.find(name)==m_descriptors.end())
		return nullptr;
	else
		return m_descriptors[name];

}

shared_ptr<CAbstractDescriptor> CFeature::GetDescriptor(size_t no) {

	if(no>=m_descriptors.size())
		return nullptr;

	map<string,shared_ptr<CAbstractDescriptor> >::iterator it = m_descriptors.begin();

	for(size_t i=0; i<no; i++)
		it++;

	return it->second;

}

string CFeature::GetDescriptorName(size_t no) {

    if(no>=m_descriptors.size())
        return nullptr;

    map<string,shared_ptr<CAbstractDescriptor> >::iterator it = m_descriptors.begin();

    for(size_t i=0; i<no; i++)
        it++;

    return it->first;

}




bool CFeature::HasDescriptor(const char* name) {

	if(!m_descriptors.empty())
		return m_descriptors.find(name)!=m_descriptors.end();
	else
		return false;

}

vec CFeature::GetLocationAtNativeScale() {

	vec result = m_location;

	result.Scale(pow(2,m_scale));

	return result;

}

bool CFeature::SaveToFile(const char* filename, list<CFeature>& features) {

    list<CFeature>::iterator it;

    ofstream out(filename);

    if(!out.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    out << "# created by r4r_motion" << endl;

    // write number of features
    out << features.size() << endl;

    for(it=features.begin(); it!=features.end(); it++)
        out << *it;

    out.close();

}

bool CFeature::OpenFromFile(const char* filename, std::list<CFeature>& features, int type) {

    ifstream in(filename);

    if(!in.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    // read comment line
    string comment;
    getline(in,comment);

    // read number of features
    size_t n;
    in >> n;
    in.get();

    for(size_t i=0; i<n; i++) {

        CFeature x;

        // read location
        in >> x.m_location(0);
        in >> x.m_location(1);
        in.get();

        // scale
        in >> x.m_scale;
        in.get();

        // read quality
        in >> x.m_quality;
        in.get();

        // read number of descriptors
        size_t nod;
        in >> nod;
        in.get();

        for(size_t j=0; j<nod; j++) {

            // name
            string name;
            getline(in,name);

            // read all as floats, FIXME: add other formats
            switch(type) {

            case F4S:
            {

                CDescriptor<vecf>* pdesc = new CDescriptor<vecf>();
                pdesc->Read(in);
                x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

                break;

            }
            }

        }

        features.push_back(x);

    }

    in.close();

}

bool CFeature::SaveDescriptors(const char* filename, list<CFeature>& features, const char* name, const char* comment, int type, size_t t0) {

    string fndata;

    // check if file exists
    ifstream in(filename);
    bool exists = in.good();

    if(exists) {

        string comment;
        getline(in,comment);
        getline(in,fndata);

    }

    in.close();

    // open header file
    ofstream headers(filename,ios_base::app);

    if(!headers.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    // if the file is empty, create header
    if(!exists) {

        headers << "# created by r4r_motion" << endl;

        stringstream ss;
        ss << filename << ".dat";
        fndata = ss.str();
        headers << fndata;

    }

    list<CFeature>::iterator it;

    for(it = features.begin(); it!=features.end(); it++) {

        // check if descriptor exists
        if(it->HasDescriptor(name)) {

            // break line first
            headers << endl;

            // create header and write it
            CDescriptorFileHeader header(*it);
            header.SetDescriptor(*it,name,type);
            header.SetTime(t0);
            header.SetComment(comment);

            headers << header;
        }


    }

    headers.close();

    // now write data
    ofstream data(fndata,ios_base::app);

    if(!data.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    for(it = features.begin(); it!=features.end(); it++) {


        if(it->HasDescriptor(name)) {

            CAbstractDescriptor* padesc = it->GetDescriptor(name).get();
            CDescriptor<vecf>* pdesc = static_cast<CDescriptor<vecf>*>(padesc);              // make distinction between data types
            vecf& container = pdesc->Get();
            data.write((char*)(container.Data()),sizeof(float)*container.NElems());

        }

    }

    data.close();

    return 0;

}

float* CFeature::LoadDescriptors(const char* filename, vector<CDescriptorFileHeader>& headers) {

    string fndata;

    // check if file exists
    ifstream in(filename);

    if(!in.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return NULL;

    }

    if(in.good()) {

        string comment;
        getline(in,comment);
        getline(in,fndata);

    }

    // read headers
    size_t nelems, counter;
    counter = 0;
    nelems = 0;

    while(in.good()) {

        CDescriptorFileHeader header;
        in >> header;
        in.get();
        headers.push_back(header);

        nelems += header.NElems();                      // we assume that the containers are floats
        counter ++;

    }

    in.close();

    ifstream ind(fndata.c_str());

    if(!ind.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return NULL;

    }

    // create buffer or memory map
    float* data = new float[nelems];

    // read into buffer
    ind.read((char*)data,sizeof(float)*nelems);

    ind.close();

    cout << "Read " << counter << " descriptors." << endl;

    return data;

}


}

