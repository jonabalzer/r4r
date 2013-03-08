#include "dagg.h"
#include "types.h"


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
CFeature CDescriptorAggregator<Array>::CopyFeature(CFeature x) {

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

template class CDescriptorAggregator<matf>;
template class CSubsampleAggregator<matf>;
template class CInitFrameAggregator<matf>;

}
