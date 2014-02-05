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
#include "rbuffer.h"

#include <stdio.h>
#include <assert.h>

using namespace cv;
using namespace std;

namespace R4R {

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint():
    m_location(),
    m_scale(0),
    m_quality(0),
    m_descriptors() {}

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint(std::initializer_list<T> ilist):
    m_location(),
    m_scale(0),
    m_quality(0),
    m_descriptors() {

    u_int i = 0;

    for(auto it = ilist.begin(); it!=ilist.end(); ++i, ++it)
        m_location(i) = *it;

}

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint(const CVector<T,n>& location):
    m_location(location),
    m_scale(0),
    m_quality(0),
    m_descriptors() {}

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint(const CVector<T,n>& location, float scale, T quality):
    m_location(location),
    m_scale(scale),
    m_quality(quality),
    m_descriptors() {}

template<typename T,u_int n>
CInterestPoint<T,n>::CInterestPoint(const CInterestPoint<T,n>& x, string name):
    m_location(x.m_location),
    m_scale(x.m_scale),
    m_quality(x.m_quality),
    m_descriptors() {

    unordered_map<string,shared_ptr<CAbstractDescriptor> >::const_iterator it;

    for(it=x.m_descriptors.begin(); it!=x.m_descriptors.end(); ++it) {

        if(it->first==name) {

            shared_ptr<CAbstractDescriptor> pdesc = x.GetDescriptor(name.c_str());
            this->AttachDescriptor(name.c_str(),pdesc);

        }

    }

}




template<typename T,u_int n>
CInterestPoint<T,n>::~CInterestPoint() {

    unordered_map<string,shared_ptr<CAbstractDescriptor> >::iterator it;

    for(it=m_descriptors.begin(); it!=m_descriptors.end(); ++it)
        it->second.reset();

}

template<typename T,u_int n>
bool CInterestPoint<T,n>::operator==(const CInterestPoint<T,n>& x) {

    return m_location==x.m_location && m_scale==x.m_scale;

}

template<typename T,u_int n>
void CInterestPoint<T,n>::AttachDescriptor(const char* id, shared_ptr<CAbstractDescriptor> descriptor) {

    m_descriptors.insert(pair<string,shared_ptr<CAbstractDescriptor> >(id,descriptor));

}

template<typename T,u_int n>
bool CInterestPoint<T,n>::HasDescriptor(const char* name) const {

    if(!m_descriptors.empty())
        return m_descriptors.find(name)!=m_descriptors.end();
    else
        return false;

}

template<typename T,u_int n>
CVector<T,n> CInterestPoint<T,n>::GetLocationAtNativeScale() const {

    return pow(2,m_scale)*m_location;

}

template<typename T,u_int n>
const shared_ptr<CAbstractDescriptor>& CInterestPoint<T, n>::GetDescriptor(const char* name) const {

    if(m_descriptors.find(name)==m_descriptors.end())
        return nullptr;
    else
        return m_descriptors.at(name);

}

template<typename T,u_int n>
const shared_ptr<CAbstractDescriptor>& CInterestPoint<T,n>::GetDescriptor(u_int no) const {

    if(no>=m_descriptors.size())
        return nullptr;

    unordered_map<string,shared_ptr<CAbstractDescriptor> >::const_iterator it = m_descriptors.begin();

    for(size_t i=0; i<no; i++)
        it++;

    return it->second;

}

template<typename T,u_int n>
string CInterestPoint<T,n>::GetDescriptorName(u_int no) {

    if(no>=m_descriptors.size())
        return nullptr;

    unordered_map<string,shared_ptr<CAbstractDescriptor> >::iterator it = m_descriptors.begin();

    for(size_t i=0; i<no; i++)
        it++;

    return it->first;

}

template<typename U,u_int m>
ostream& operator<<(ostream& os, const CInterestPoint<U,m>& x) {

    os << "Location: " << x.m_location << endl;
    os << "Scale: " << x.m_scale << endl;
    os << "Quality: " << x.m_quality << endl;
    os << "Number of descriptors: " << x.NoDescriptors();

    return os;

}

template ostream& operator<<(ostream& os, const CInterestPoint<float,2>& x);
template ostream& operator<<(ostream& os, const CInterestPoint<float,3>& x);

template<typename U,u_int m>
ofstream& operator<<(ofstream& os, const CInterestPoint<U,m>& x) {

    // location, scale, and quality
    for(uint i=0; i<m; i++) {

        os << x.m_location.Get(i);

        if(i<m-1)
            os << " ";
        else
            os << endl;

    }

    os << x.m_scale << endl;
    os << x.m_quality << endl;

    // number of descriptors that follow
    u_int nod = x.NoDescriptors();
    if(nod>0)
        os << x.NoDescriptors() << endl;
    else
        os << x.NoDescriptors();

    unordered_map<string,shared_ptr<CAbstractDescriptor> >::const_iterator it;

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

template ofstream& operator<<(ofstream& os, const CInterestPoint<float,2>& x);
template ofstream& operator<<(ofstream& os, const CInterestPoint<float,3>& x);

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

template ifstream& operator>>(ifstream& is, CInterestPoint<float,2>& x);

template<typename T,u_int n>
template<template<class U, class Allocator = std::allocator<U> > class Container>
bool CInterestPoint<T,n>::SaveToFile(const char* filename, const Container<CInterestPoint<T,n> >& features, const char* comment) {

    typename Container<CInterestPoint<T,n> >::const_iterator it;

    ofstream out(filename);

    if(!out.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    if(comment==nullptr)
        out << "# created by r4r_motion" << endl;
    else
        out << "# " << comment << endl;

    // write number of features (for matlab import)
    out << features.size() << endl;

    size_t counter = 0;

    for(it=features.begin(); it!=features.end(); ++it, ++counter) {

        // write feature itself
        out << *it;

        // do not break line after the last feature
        if(counter<features.size()-1)
            out << endl;

    }

    out.close();

    return 0;

}

template bool CInterestPoint<float,2>::SaveToFile(const char* filename, const vector<CInterestPoint<float,2> >& features, const char* comment);
template bool CInterestPoint<float,2>::SaveToFile(const char* filename, const list<CInterestPoint<float,2> >& features, const char* comment);
template bool CInterestPoint<float,2>::SaveToFile(const char* filename, const CRingBuffer<CInterestPoint<float,2> >& features, const char* comment);
template bool CInterestPoint<float,3>::SaveToFile(const char* filename, const vector<CInterestPoint<float,3> >& features, const char* comment);
template bool CInterestPoint<float,3>::SaveToFile(const char* filename, const list<CInterestPoint<float,3> >& features, const char* comment);
template bool CInterestPoint<float,3>::SaveToFile(const char* filename, const CRingBuffer<CInterestPoint<float,3> >& features, const char* comment);


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

    int nis;
    in >> nis;
    in.get();

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
template class CInterestPoint<float,3>;

}

