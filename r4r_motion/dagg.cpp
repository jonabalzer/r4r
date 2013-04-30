#include "dagg.h"
#include "types.h"
#include "descriptor.h"
#include "factor.h"
#include "kernels.h"

namespace R4R {

template <class Array>
CDescriptorAggregator<Array>::CDescriptorAggregator(CTracker* tracker, const char* name):
    m_tracker(tracker),
    m_name(string(name)),
    m_aggregate() {

}

template <class Array>
void CDescriptorAggregator<Array>::Aggregate() {

    cout << "Aggregating descriptors..." << endl;

    list<shared_ptr<CTracklet> >::iterator it;

    for(size_t s=0; s<m_tracker->size(); s++) {

        for(it=m_tracker->at(s).begin(); it!=m_tracker->at(s).end(); it++) {

            AggregateTracklet((*it).get());

        }

    }
}

template <class Array>
void CDescriptorAggregator<Array>::AggregateTracklet(CTracklet* tracklet) {

    list<CFeature>::iterator it;

    for(it=tracklet->begin(); it!=tracklet->end(); it++) {

        if(it->HasDescriptor(m_name.c_str())) {

            CFeature x = CDescriptorAggregator<Array>::CopyFeature(*it);
            m_aggregate.push_back(x);

        }

    }

}

template <class Array>
bool CDescriptorAggregator<Array>::Aggregate(const char* filename, const char* comment) {

    string fndata;

    // check if file exists
    ifstream in(filename);
    bool exists = in.good();

    if(exists) {

        string usercomment;
        getline(in,usercomment);
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

    // write headers
    list<shared_ptr<CTracklet> >::iterator it;

    // go through tracklets
    for(size_t s=0; s<m_tracker->size(); s++) {

        for(it=m_tracker->at(s).begin(); it!=m_tracker->at(s).end(); it++) {

            string hash = (*it)->GetHash();

            list<CFeature>::iterator itf;

            for(itf=(*it)->begin(); itf!=(*it)->end(); itf++) {

                if(itf->HasDescriptor(m_name.c_str())) {

                    // break line first
                    headers << endl;

                    // create header and write it
                    CDescriptorFileHeader header(*itf);
                    header.SetDescriptor(*itf,m_name.c_str());
                    header.SetComment((string(comment)+string("_")+hash).c_str());

                    headers << header;

                }

            }

        }

    }

    headers.close();

    // now write data
    ofstream data(fndata,ios_base::app);

    if(!data.is_open()) {

        cout << "ERROR: Could not open file." << endl;
        return 1;

    }

    // go through tracklets again
    for(size_t s=0; s<m_tracker->size(); s++) {

        for(it=m_tracker->at(s).begin(); it!=m_tracker->at(s).end(); it++) {

            list<CFeature>::iterator itf;

            for(itf=(*it)->begin(); itf!=(*it)->end(); itf++) {

                if(itf->HasDescriptor(m_name.c_str())) {

                    shared_ptr<CAbstractDescriptor> pdesc = itf->GetDescriptor(m_name.c_str());

                    // write type
                    ETYPE type = pdesc->GetType();

                    // access to data
                    void* pdata = pdesc->GetData();

                    // write data depending on number of bytes
                    if(type==B1U || type == C1U || type == C1S)
                        data.write((char*)(pdata),sizeof(char)*pdesc->NElems());
                    else if(type==S2U || type == S2S)
                        data.write((char*)(pdata),sizeof(short)*pdesc->NElems());
                    else if(type==I4S || type==I4U || type==F4S)
                        data.write((char*)(pdata),sizeof(int)*pdesc->NElems());
                    else if(type==L8S || type== L8U || type==D8S)
                        data.write((char*)(pdata),sizeof(double)*pdesc->NElems());

                }

            }

        }

    }

    data.close();

    return 0;

}



template <class Array>
CFeature CDescriptorAggregator<Array>::CopyFeature(CFeature x) {

    // don't copy all the old descriptors
    CFeature result(x.GetLocation(),x.GetScale(),x.GetQuality());

    if(x.HasDescriptor(m_name.c_str())) {

        // get from x
        shared_ptr<CAbstractDescriptor> pdesc = x.GetDescriptor(m_name.c_str());

        // attach to result
        result.AttachDescriptor(m_name.c_str(),pdesc);

    }

    return result;

}

template <class Array>
void CInitFrameAggregator<Array>::AggregateTracklet(CTracklet* tracklet) {

    CFeature x = CDescriptorAggregator<Array>::CopyFeature(*tracklet->begin());

    if(x.HasDescriptor(m_name.c_str()))
        m_aggregate.push_back(x);

}

template <class Array>
CSubsampleAggregator<Array>::CSubsampleAggregator(CTracker* tracker, const char* name, size_t n):
    CDescriptorAggregator<Array>::CDescriptorAggregator(tracker,name),
    m_n(n) {}


template <class Array>
void CSubsampleAggregator<Array>::AggregateTracklet(CTracklet* tracklet) {

    list<CFeature>::iterator it;

    size_t counter = 0;
    for(it=tracklet->begin(); it!=tracklet->end(); it++) {

        if(counter%m_n==0 && it->HasDescriptor(m_name.c_str())) {

            CFeature x = CDescriptorAggregator<Array>::CopyFeature(*it);
            m_aggregate.push_back(x);

        }

        counter++;

    }

}


template <class Array>
CSplineInterpolationAggregator<Array>::CSplineInterpolationAggregator(CTracker* tracker, const char* name, size_t n, size_t p):
    CDescriptorAggregator<Array>::CDescriptorAggregator(tracker,name),
    m_n(n),
    m_p(p) {

}

template <class Array>
void CSplineInterpolationAggregator<Array>::AggregateTracklet(CTracklet* tracklet) {

    // access the first descriptor
    list<CFeature>::iterator it = tracklet->begin();

    shared_ptr<CDescriptor<Array> > pdesc;
    if(it->HasDescriptor(m_name.c_str()))
        pdesc = static_pointer_cast<CDescriptor<Array> >(it->GetDescriptor(m_name.c_str()));

    // get dimension of descriptors
    size_t d = pdesc->NElems();

    // create spline object
    CSplineCurve<float> curve(d,m_p,m_n);
    curve.MakeClampedUniformKnotVector(0,1);

    // set up interpolation matrix
    matf A(tracklet->size(),m_n);
    float dt = 1.0/((float)tracklet->size()-1.0);
    size_t order = m_p + 1;
    float* N = new float[order*order];
    float* knot = curve.GetKnotVector().Data().get();

    for(size_t i=0; i<tracklet->size(); i++) {

        int span = curve.GetSpan(i*dt);

        CSplineCurve<float>::EvaluateNurbsBasis(order,knot+span,i*dt,N);

        for(size_t j=0; j<order; j++)
            A(i,j) = N[j];

    }

    // compute Moore-Penrose inverse by svd, TODO: do sparse iterative SVD
    matf U, S, Vt;
    CMatrixFactorization<float>::SVD(A,U,S,Vt);

    // compute inverse of S
    size_t rank = 0;

    for(size_t i=0; i<min(S.NRows(),S.NCols()); i++) {

        if(S.Get(i,i)) {
            S(i,i) = 1.0/S.Get(i,i);
            rank++;
        }
        else
            break; // the singular values are order (number should be known)

    }

    // buffer for coefficients
    vecf tube(tracklet->size());

    // access to cv
    matf& cv = curve.GetCVData();

    // create kernel for fast inversion
    CMercerKernel<float> kernel(tracklet->size());

    // go through rows of cv (one element of the descriptor)
    for(size_t i=0; i<pdesc->NRows(); i++) {

        for(size_t j=0; j<pdesc->NCols(); j++) {

            size_t counter = 0;

            // assemble tube
            for(it=tracklet->begin(); it!=tracklet->end(); it++) {

                // downcast pointer
                pdesc = static_pointer_cast<CDescriptor<Array> >(it->GetDescriptor(m_name.c_str()));

                // get access to container, cast to float
                tube(counter) = (float)pdesc->Get().Get(i,j);

            }

            // access to solution
            vecf result = cv.GetColumn(i*pdesc->NRows()+j);

            // solve linear least-squares problem
            for(size_t k=0; k<rank; k++) {

                // project onto range of A and invert
                vecf colu = U.GetColumn(k);
                float coeff = S(k,k)*kernel.Evaluate(tube.Data().get(),colu.Data().get());

                // express in basis of range(At)
                vecf colv = Vt.GetColumn(k).Clone();
                colv.Scale(coeff);
                result.Add(colv);

            }

        }

    }

    // create new feature/descriptor pair
    CFeature x0 = *tracklet->begin();
    CFeature result(x0.GetLocation(),x0.GetScale(),x0.GetQuality());
    shared_ptr<CDescriptor<matf> > desc = shared_ptr<CDescriptor<matf> >(new CDescriptor<matf>(cv));
    result.AttachDescriptor("SPLINE",desc);
    m_aggregate.push_back(result);

}


template class CDescriptorAggregator<matf>;
template class CSplineInterpolationAggregator<matf>;
template class CSubsampleAggregator<matf>;
template class CInitFrameAggregator<matf>;
template class CDescriptorAggregator<vecf>;
template class CSubsampleAggregator<vecf>;
template class CInitFrameAggregator<vecf>;

}
