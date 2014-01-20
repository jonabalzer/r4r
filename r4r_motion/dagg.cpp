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

#include "dagg.h"
#include "types.h"
#include "descriptor.h"
#include "factor.h"
#include "kernels.h"

namespace R4R {

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
CDescriptorAggregator<Array,Container>::CDescriptorAggregator(CTracker<Container>* tracker, const char* name):
    m_tracker(tracker),
    m_name(string(name)),
    m_aggregate() {

}

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
void CDescriptorAggregator<Array,Container>::Aggregate() {

    cout << "Aggregating descriptors..." << endl;

    const Container<shared_ptr<CTracklet> >& data = m_tracker->GetData();
    typename Container<shared_ptr<CTracklet> >::const_iterator it;

    for(it=data.begin(); it!=data.end(); it++) {

        cout << ".";

        AggregateTracklet((*it).get());

    }

    cout << endl;

}

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
void CDescriptorAggregator<Array,Container>::AggregateTracklet(CTracklet* tracklet) {

    // first compress and reverse order
    tracklet->CompressAndReverse();

    // get access to the data
    const vector<imfeature>& data = tracklet->GetData();
    vector<imfeature>::const_iterator it = data.begin();

    // keep the first feature as reference
    imfeature x0 = *it;

    for(it; it!= data.end(); it++) {

        if(it->HasDescriptor(m_name.c_str())) {

            imfeature x = CDescriptorAggregator<Array,Container>::CopyFeature(*it);

            // copy properties of reference feature
            x.m_scale = x0.m_scale;
            x.m_quality = x0.m_quality;
            x.m_location(0) = x0.m_location.Get(0);
            x.m_location(1) = x0.m_location.Get(1);

            m_aggregate.push_back(x);

        }

    }

}

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
imfeature CDescriptorAggregator<Array,Container>::CopyFeature(const imfeature& x) const {

    // don't copy all the old descriptors
    const vec2f& loc = x.GetLocation();
    imfeature result(loc,x.GetScale(),x.GetQuality());

    if(x.HasDescriptor(m_name.c_str())) {

        // get from x
        shared_ptr<CAbstractDescriptor> pdesc = x.GetDescriptor(m_name.c_str());

        // attach to result
        result.AttachDescriptor(m_name.c_str(),pdesc);

    }

    return result;

}

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
void CInitFrameAggregator<Array,Container>::AggregateTracklet(CTracklet* tracklet) {

    // first compress and reverse order
    tracklet->CompressAndReverse();

    // get access to the data
    const vector<imfeature>& data = tracklet->GetData();
    vector<imfeature>::const_iterator it = data.begin();

    imfeature x = CDescriptorAggregator<Array,Container>::CopyFeature(*it);

    if(x.HasDescriptor(m_name.c_str()))
        m_aggregate.push_back(x);

}

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
CSubsampleAggregator<Array,Container>::CSubsampleAggregator(CTracker<Container>* tracker, const char* name, size_t n):
    CDescriptorAggregator<Array,Container>::CDescriptorAggregator(tracker,name),
    m_n(n) {}


template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
void CSubsampleAggregator<Array,Container>::AggregateTracklet(CTracklet* tracklet) {

    // first compress and reverse order
    tracklet->CompressAndReverse();

    // get access to the data
    const vector<imfeature>& data = tracklet->GetData();
    vector<imfeature>::const_iterator it = data.begin();

    // keep the first feature as reference
    imfeature x0 = *it;
    size_t counter = 0;

    for(it; it!= data.end(); it++) {

        if(counter%m_n==0 && it->HasDescriptor(m_name.c_str())) {

            imfeature x = CDescriptorAggregator<Array,Container>::CopyFeature(*it);

            // copy properties of reference feature
            x.m_scale = x0.m_scale;
            x.m_quality = x0.m_quality;
            x.m_location(0) = x0.m_location.Get(0);
            x.m_location(1) = x0.m_location.Get(1);

            m_aggregate.push_back(x);

        }

        counter++;

    }

}


template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
void CMeanAggregator<Array,Container>::AggregateTracklet(CTracklet* tracklet) {

    // first compress and reverse order
    tracklet->CompressAndReverse();

    // get access to the data
    const vector<imfeature>& data = tracklet->GetData();
    vector<imfeature>::const_iterator it = data.begin();

    shared_ptr<CDescriptor<Array> > pdesc;
    if(it->HasDescriptor(m_name.c_str()))
        pdesc = static_pointer_cast<CDescriptor<Array> >(it->GetDescriptor(m_name.c_str()));
    else {

        cerr << "ERROR: Descriptor " << m_name << " not found!" << endl;
        return;

    }

    Array mean(pdesc->NRows(),pdesc->NCols());

    size_t counter = 0;

    for(it; it!=data.end(); it++, counter++) {

        // only add and count features that have the descripor
        if(it->HasDescriptor(m_name.c_str())) {

            pdesc = static_pointer_cast<CDescriptor<Array> >(it->GetDescriptor(m_name.c_str()));
            mean = mean + pdesc->Get();

        }

    }

    if(counter>0)
        mean.Scale(1.0/(double)counter);

    // create new feature/descriptor pair
    imfeature x0 = *data.begin();
    const vec2f& loc = x0.GetLocation();
    imfeature result(loc,x0.GetScale(),x0.GetQuality());
    shared_ptr<CDescriptor<matf> > desc = shared_ptr<CDescriptor<matf> >(new CDescriptor<matf>(mean));
    result.AttachDescriptor(m_name.c_str(),desc);
    m_aggregate.push_back(result);

}


template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
CSplineInterpolationAggregator<Array,Container>::CSplineInterpolationAggregator(CTracker<Container>* tracker, const char* name, size_t n, size_t p):
    CDescriptorAggregator<Array,Container>::CDescriptorAggregator(tracker,name),
    m_n(n),
    m_p(p) {

}

template <class Array,template<class T,class Allocator = std::allocator<T> > class Container>
void CSplineInterpolationAggregator<Array,Container>::AggregateTracklet(CTracklet* tracklet) {

    // first compress and reverse order
    tracklet->CompressAndReverse();

    // get access to the data
    const vector<imfeature>& data = tracklet->GetData();

    // if the tracklet is too short disregard it, A must have full rank
    if(data.size()<m_n)
        return;

    // access the first descriptor
    vector<imfeature>::const_iterator it = data.begin();

    shared_ptr<CDescriptor<Array> > pdesc;
    if(it->HasDescriptor(m_name.c_str()))
        pdesc = static_pointer_cast<CDescriptor<Array> >(it->GetDescriptor(m_name.c_str()));
    else
        return;

    // get dimension of descriptors
    size_t d = pdesc->NElems();

    // create spline object
    CSplineCurve<float> curve(d,m_p,m_n);
    curve.MakeClampedUniformKnotVector(0,1);

    // set up interpolation matrix
    matf A(data.size(),m_n);
    float dt = 1.0/((float)data.size()-1.0);
    size_t order = m_p + 1;
    float* N = new float[order*order];
    float* knot = curve.GetKnotVector().Data().get();

    for(size_t i=0; i<data.size(); i++) {

        int span = curve.GetSpan(i*dt);

        CSplineCurve<float>::EvaluateNurbsBasis(order,knot+span,i*dt,N);

        for(size_t j=0; j<order; j++)
            A(i,span+j) = N[j];

    }

    delete [] N;

    // compute Moore-Penrose inverse by svd
    matf U(A.NRows(),A.NRows());
    vecf s(min(A.NRows(),A.NCols()));
    matf Vt(A.NCols(),A.NCols());

    CMatrixFactorization<float>::SVD(A,U,s,Vt);

    // compute inverse of s
    size_t rank = 0;

    for(size_t i=0; i<s.NElems(); i++) {

        if(s.Get(i)>0) {
            s(i) = 1.0/s.Get(i);
            rank++;
        }
        else
            break; // the singular values are ordered (number should be known)

    }

    // buffer for coefficients
    vecf tube(data.size());

    // access to cv
    matf& cv = curve.GetCVData();

    // create kernel for fast inversion
    CMercerKernel<float> kernel(data.size());

    // iterate through all pixels
    for(size_t i=0; i<pdesc->NRows(); i++) {

        for(size_t j=0; j<pdesc->NCols(); j++) {

            size_t counter = 0;

            // assemble tube for the pixel i,j
            for(it=data.begin(); it!=data.end(); it++) {

                // downcast pointer
                pdesc = static_pointer_cast<CDescriptor<Array> >(it->GetDescriptor(m_name.c_str()));

                // get access to container, cast to float
                tube(counter) = (float)pdesc->Get().Get(i,j);

                counter++;

            }

            // access to solution
            vecf result(1,m_n);

            // solve linear least-squares problem
            for(size_t k=0; k<rank; k++) {

                // project onto range of A and invert
                vecf colu = U.GetColumn(k);

                float coeff = s.Get(k)*kernel.Evaluate(tube.Data().get(),colu.Data().get());

                // express in basis of range(At)
                vecf colv = Vt.GetRow(k);           // col of Vt = row of V, no need to clone here
                colv.Scale(coeff);
                result.Add(colv);

            }

            // copy result into control point matrix
            for(size_t l=0; l<result.NElems(); l++)
                cv(i*pdesc->NCols()+j,l) = result.Get(l);

        }

    }

    // create new feature/descriptor pair
    imfeature x0 = *data.begin();
    vec2f loc = x0.GetLocation();
    imfeature result(loc,x0.GetScale(),x0.GetQuality());
    shared_ptr<CDescriptor<matf> > desc = shared_ptr<CDescriptor<matf> >(new CDescriptor<matf>(cv));
    result.AttachDescriptor(m_name.c_str(),desc);
    m_aggregate.push_back(result);

}


template class CDescriptorAggregator<matf,list>;
template class CSplineInterpolationAggregator<matf,list>;
template class CMeanAggregator<matf,list>;
template class CSubsampleAggregator<matf,list>;
template class CInitFrameAggregator<matf,list>;
template class CDescriptorAggregator<vecf,list>;
template class CSubsampleAggregator<vecf,list>;
template class CInitFrameAggregator<vecf,list>;
template class CSplineInterpolationAggregator<vecf,list>;
template class CMeanAggregator<vecf,list>;

}
