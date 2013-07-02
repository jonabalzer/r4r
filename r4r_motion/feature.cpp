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

CFeature::CFeature(const CFeature& x):
    m_location(2),
    m_scale(x.m_scale),
    m_quality(x.m_quality),
    m_descriptors(x.m_descriptors) {

    m_location(0) = x.m_location.Get(0);
    m_location(1) = x.m_location.Get(1);

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
    os << "Number of descriptors: " << x.NoDescriptors();

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
    size_t nod = x.NoDescriptors();
    if(nod>0)
        os << x.NoDescriptors() << endl;
    else
        os << x.NoDescriptors();

    map<string,shared_ptr<CAbstractDescriptor> >::iterator it;

	// print all descriptors
    size_t counter = 0;
    for(it = x.m_descriptors.begin(); it!=x.m_descriptors.end(); it++, counter++) {

		os << it->first.c_str() << endl;
        it->second->Write(os);

        if(counter<nod-1)       // only break line between the descriptors
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

        // get size of descriptor
        size_t nrows, ncols;
        is >> nrows;
        is >> ncols;

        // read type, TODO: map type name to type
        int temp;
        is >> temp;
        ETYPE type = (ETYPE)temp;
        is.get();

        // treat different data types separately
        switch(type) {

        case ETYPE::B1U:
        {

            CDenseArray<bool> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(bool)*container.NElems());

            CDescriptor<CDenseArray<bool> >* pdesc = new CDescriptor<CDenseArray<bool> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::I4S:
        {

            CDenseArray<int> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(int)*container.NElems());

            CDescriptor<CDenseArray<int> >* pdesc = new CDescriptor<CDenseArray<int> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::F4S:
        {

            CDenseArray<float> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(float)*container.NElems());

            CDescriptor<CDenseArray<float> >* pdesc = new CDescriptor<CDenseArray<float> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::L8U:
        {

            CDenseArray<size_t> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(size_t)*container.NElems());

            CDescriptor<CDenseArray<size_t> >* pdesc = new CDescriptor<CDenseArray<size_t> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::D8S:
        {

            CDenseArray<double> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(double)*container.NElems());

            CDescriptor<CDenseArray<double> >* pdesc = new CDescriptor<CDenseArray<double> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }


        default:
            return is;

        }

        is.get();


    }



	return is;

}

bool CFeature::operator!=(CFeature& x) {

	vec pt0 = this->GetLocation();
	vec pt1 = x.GetLocation();

	return pt0(0)!=pt1(0) && pt0(1)!=pt1(1);

}

//CFeature::~CFeature() {

//	DeleteDescriptors();

//}

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

bool CFeature::SaveToFile(const char* filename, list<CFeature>& features, const char* comment) {

    list<CFeature>::iterator it;

    ofstream out(filename);

    if(!out.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    if(comment==nullptr)
        out << "# created by r4r_motion" << endl;
    else
        out << "# " << comment << endl;

    size_t counter = 0;

    for(it=features.begin(); it!=features.end(); it++, counter++) {
        out << *it;

        if(counter<features.size()-1)
            out << endl;

    }

    out.close();

    return 0;

}

int CFeature::OpenFromFile(const char* filename, std::vector<CFeature>& features, string& comment) {

    ifstream in(filename);

    if(!in.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return -1;

    }

    // read comment line
    //string comment;
    getline(in,comment);
    comment = comment.substr(2,comment.size()-1);

    int n = 0;

    while(true) {

        if(in.good()) {

            CFeature x;
            in >> x;

            n++;
            features.push_back(x);

        }
        else
            break;
    }

    in.close();

    return n;

}

//bool CFeature::SaveDescriptors(const char* filename, list<CFeature>& features, const char* name, const char* comment, ETYPE type, size_t t0) {

//    string fndata;

//    // check if file exists
//    ifstream in(filename);
//    bool exists = in.good();

//    if(exists) {

//        string usercomment;
//        getline(in,usercomment);
//        getline(in,fndata);

//    }

//    in.close();

//    // open header file
//    ofstream headers(filename,ios_base::app);

//    if(!headers.is_open()) {

//        cout << "ERROR: Could not open file." << endl;
//        return 1;

//    }

//    // if the file is empty, create header
//    if(!exists) {

//        headers << "# created by r4r_motion" << endl;

//        stringstream ss;
//        ss << filename << ".dat";
//        fndata = ss.str();
//        headers << fndata;

//    }

//    list<CFeature>::iterator it;

//    for(it = features.begin(); it!=features.end(); it++) {

//        // check if descriptor exists
//        if(it->HasDescriptor(name)) {

//            // break line first
//            headers << endl;

//            // create header and write it
//            CDescriptorFileHeader header(*it);      // sets feature properties
//            header.SetDescriptor(*it,name);         // sets descriptor props: size, name, and type
//            header.SetComment(comment);             // set comments, combines object description and tracklet hash

//            headers << header;

//        }

//    }

//    headers.close();

//    // now write data
//    ofstream data(fndata,ios_base::app);

//    if(!data.is_open()) {

//        cout << "ERROR: Could not open file." << endl;
//        return 1;

//    }

//    for(it = features.begin(); it!=features.end(); it++) {

//        if(it->HasDescriptor(name)) {

//            CAbstractDescriptor* padesc = it->GetDescriptor(name).get();
//            CDescriptor<vecf>* pdesc = static_cast<CDescriptor<vecf>*>(padesc);              // make distinction between data types
//            vecf& container = pdesc->Get();
//            data.write((char*)(container.Data()),sizeof(float)*container.NElems());

//        }

//    }

//    data.close();

//    return 0;

//}

//float* CFeature::LoadDescriptors(const char* filename, vector<CDescriptorFileHeader>& headers) {

//    string fndata;

//    // check if file exists
//    ifstream in(filename);

//    if(!in.is_open()) {

//        cout << "ERROR: Could not open file." << endl;
//        return NULL;

//    }

//    if(in.good()) {

//        string comment;
//        getline(in,comment);
//        getline(in,fndata);

//    }

//    // read headers
//    size_t nelems, counter;
//    counter = 0;
//    nelems = 0;

//    while(in.good()) {

//        CDescriptorFileHeader header;
//        in >> header;
//        in.get();
//        headers.push_back(header);

//        nelems += header.NElems();                      // we assume that the containers are floats
//        counter ++;

//    }

//    in.close();

//    ifstream ind(fndata.c_str());

//    if(!ind.is_open()) {

//        cout << "ERROR: Could not open file." << endl;
//        return NULL;

//    }

//    // create buffer or memory map
//    float* data = new float[nelems];

//    // read into buffer
//    ind.read((char*)data,sizeof(float)*nelems);

//    ind.close();

//    cout << "Read " << counter << " descriptors." << endl;

//    return data;

//}

//bool CFeature::LoadDescriptors(const char *filename, std::vector<CDescriptorFileHeader>& headers) { //, void *data, ETYPE &type) {

//    return 0;
//}


void* CFeature::LoadDescriptors(const char* filename, std::vector<CDescriptorFileHeader>& headers, ETYPE& type, bool preview) {

    string fndata;

    // check if file exists
    ifstream in(filename);

    if(!in.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return nullptr;

    }

    if(in.good()) {

        string comment;
        getline(in,comment);
        getline(in,fndata);

    }

    // read first header
    size_t nelems0, nelems, counter;
    ETYPE type0;

    counter = 0;
    nelems = 0;

    if(in.good()) {

        CDescriptorFileHeader header;
        in >> header;
        in.get();
        headers.push_back(header);

        nelems = header.NElems();
        nelems0 = nelems;
        type0 = header.m_type;

        counter ++;

    }

    while(in.good()) {

        CDescriptorFileHeader header;
        in >> header;
        in.get();

        if(header.NElems()!=nelems0 || header.m_type != type0) {

            cout << "ERROR: Inconsitent size or type." << endl;
            return nullptr;

        }

        headers.push_back(header);
        nelems += header.NElems();
        counter ++;

    }

    in.close();

    if(preview)
        return nullptr;

    ifstream ind(fndata.c_str());

    if(!ind.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return nullptr;

    }

    // create buffer or memory map
    void* data;

    switch(type0) {

    case ETYPE::B1U:
    {

        type = ETYPE::B1U;

        data = new bool[nelems];

        // read into buffer
        ind.read((char*)data,sizeof(bool)*nelems);

        break;

    }

    case ETYPE::I4S:
    {

        type = ETYPE::I4S;

        data = new int[nelems];

        ind.read((char*)data,sizeof(int)*nelems);

        break;

    }


    case ETYPE::F4S:
    {

        type = ETYPE::F4S;

        data = new float[nelems];

        ind.read((char*)data,sizeof(float)*nelems);

        break;

    }

    case ETYPE::L8U:
    {

        type = ETYPE::L8U;

        data = new size_t[nelems];

        ind.read((char*)data,sizeof(size_t)*nelems);

        break;

    }

    case ETYPE::D8S:
    {

        type = ETYPE::D8S;

        data = new double[nelems];

        ind.read((char*)data,sizeof(double)*nelems);

        break;

    }

    default:
        break;


    }

    ind.close();

    cout << "Read " << counter << " descriptors." << endl;

    return data;

}

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint():
    m_location(),
    m_scale(0),
    m_quality(0) {}

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint(CVector<T,n>& location, u_int scale, T quality):
    m_location(location),
    m_scale(scale),
    m_quality(quality) {}


template<typename T,u_int n>
bool CInterestPoint<T,n>::operator==(const CInterestPoint<T,n>& x) {

    return m_location==x.m_location && m_scale==x.m_scale;

}

template<typename T,u_int n>
void CInterestPoint<T,n>::AttachDescriptor(const char* id, shared_ptr<CAbstractDescriptor> descriptor) {

    m_descriptors.insert(pair<string,shared_ptr<CAbstractDescriptor> >(id,descriptor));

}

template<typename T,u_int n>
bool CInterestPoint<T,n>::HasDescriptor(const char* name) {

    if(!m_descriptors.empty())
        return m_descriptors.find(name)!=m_descriptors.end();
    else
        return false;

}

template<typename T,u_int n>
CVector<T,n> CInterestPoint<T,n>::GetLocationAtNativeScale() {

    return pow(2,m_scale)*m_location;

}

template<typename T,u_int n>
shared_ptr<CAbstractDescriptor> CInterestPoint<T,n>::GetDescriptor(const char* name) {

    if(m_descriptors.find(name)==m_descriptors.end())
        return nullptr;
    else
        return m_descriptors[name];

}

template<typename T,u_int n>
shared_ptr<CAbstractDescriptor> CInterestPoint<T,n>::GetDescriptor(u_int no) {

    if(no>=m_descriptors.size())
        return nullptr;

    map<string,shared_ptr<CAbstractDescriptor> >::iterator it = m_descriptors.begin();

    for(size_t i=0; i<no; i++)
        it++;

    return it->second;

}

template<typename T,u_int n>
string CInterestPoint<T,n>::GetDescriptorName(u_int no) {

    if(no>=m_descriptors.size())
        return nullptr;

    map<string,shared_ptr<CAbstractDescriptor> >::iterator it = m_descriptors.begin();

    for(size_t i=0; i<no; i++)
        it++;

    return it->first;

}

template<typename U,u_int m>
ostream& operator<<(ostream& os, CInterestPoint<U,m>& x) {

    os << "Location: " << x.m_location << endl;
    os << "Scale: " << x.m_scale << endl;
    os << "Quality: " << x.m_quality << endl;
    os << "Number of descriptors: " << x.NoDescriptors();

    return os;

}

template<typename U,u_int m>
ofstream& operator<<(ofstream& os, CInterestPoint<U,m>& x) {

    // location, scale, and quality
    os << x.m_location << endl;
    os << x.m_scale << endl;
    os << x.m_quality << endl;

    // number of descriptors that follow
    u_int nod = x.NoDescriptors();
    if(nod>0)
        os << x.NoDescriptors() << endl;
    else
        os << x.NoDescriptors();

    map<string,shared_ptr<CAbstractDescriptor> >::iterator it;

    // write all descriptors
    u_int counter = 0;
    for(it = x.m_descriptors.begin(); it!=x.m_descriptors.end(); it++, counter++) {

        os << it->first.c_str() << endl;
        it->second->Write(os);

        if(counter<nod-1)       // only break line between the descriptors
            os << endl;

    }

    return os;

}

template<typename U,u_int m>
ifstream& operator>>(ifstream& is, CInterestPoint<U,m>& x) {

    // read location
    is >> x.m_location;
    is.get();

    // scale
    is >> x.m_scale;
    is.get();

    // read quality
    is >> x.m_quality;
    is.get();

    // read number of descriptors
    u_int nod;
    is >> nod;
    is.get();

    for(size_t k=0; k<nod; k++) {

        // name
        string name;
        getline(is,name);

        // get size of descriptor
        size_t nrows, ncols;
        is >> nrows;
        is >> ncols;

        // read type
        int temp;
        is >> temp;
        ETYPE type = (ETYPE)temp;
        is.get();

        // treat different data types separately
        switch(type) {

        case ETYPE::B1U:
        {

            CDenseArray<bool> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(bool)*container.NElems());

            CDescriptor<CDenseArray<bool> >* pdesc = new CDescriptor<CDenseArray<bool> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::I4S:
        {

            CDenseArray<int> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(int)*container.NElems());

            CDescriptor<CDenseArray<int> >* pdesc = new CDescriptor<CDenseArray<int> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::F4S:
        {

            CDenseArray<float> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(float)*container.NElems());

            CDescriptor<CDenseArray<float> >* pdesc = new CDescriptor<CDenseArray<float> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::L8U:
        {

            CDenseArray<size_t> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(size_t)*container.NElems());

            CDescriptor<CDenseArray<size_t> >* pdesc = new CDescriptor<CDenseArray<size_t> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        case ETYPE::D8S:
        {

            CDenseArray<double> container(nrows,ncols);
            is.read((char*)container.Data().get(),sizeof(double)*container.NElems());

            CDescriptor<CDenseArray<double> >* pdesc = new CDescriptor<CDenseArray<double> >(container);
            x.AttachDescriptor(name.c_str(),shared_ptr<CAbstractDescriptor>(pdesc));

            break;

        }

        default:
            return is;

        }

        is.get();

    }

    return is;

}

template<typename T,u_int n>
bool CInterestPoint<T,n>::SaveToFile(const char* filename, std::list<CInterestPoint<T,n> >& features, const char* comment) {

    typename list<CInterestPoint<T,n> >::iterator it;

    ofstream out(filename);

    if(!out.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    if(comment==nullptr)
        out << "# created by r4r_motion" << endl;
    else
        out << "# " << comment << endl;

    size_t counter = 0;

    for(it=features.begin(); it!=features.end(); it++, counter++) {

        // write feature itself
        out << *it;

        // do not break line after the last feature
        if(counter<features.size()-1)
            out << endl;

    }

    out.close();

    return 0;

}

template<typename T,u_int n>
int CInterestPoint<T,n>::LoadFromFile(const char* filename, std::vector<CInterestPoint<T,n> >& features, string& comment) {

    ifstream in(filename);

    if(!in.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return -1;

    }

    // read comment line
    getline(in,comment);
    comment = comment.substr(2,comment.size()-1);

    int no = 0;

    while(true) {

        if(in.good()) {

            CInterestPoint<T,n> x;
            in >> x;

            no++;
            features.push_back(x);

        }
        else
            break;
    }

    in.close();

    return no;

}


template class CInterestPoint<float,2>;
template ostream& operator<<(ostream& os, CInterestPoint<float,2>& x);
template ofstream& operator<<(ofstream& os, CInterestPoint<float,2>& x);
template ifstream& operator>>(ifstream& is, CInterestPoint<float,2>& x);

}

