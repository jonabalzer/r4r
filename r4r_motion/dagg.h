#ifndef DAGG_H
#define DAGG_H

#include "tracker.h"
#include <list>

namespace R4R {

template <class Array>
class CDescriptorAggregator {

public:

    //! Constructor.
    CDescriptorAggregator(CTracker* tracker, const char* name);

    //! Aggregates over the entire tracker.
    void Aggregate();

    //! Aggregates over the entire tracker into a file.
    virtual bool Aggregate(const char* filename, const char* comment = nullptr);

    //! Access to the aggregates.
    list<CFeature>& Get() { return m_aggregate; };

protected:

    //! Aggregates a single tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    //! Copies the feature including only the descriptors specified by #m_name.
    CFeature CopyFeature(CFeature x);

    CTracker* m_tracker;                //! tracker to aggregate over
    string m_name;                      //! name of the descriptor to aggregate
    list<CFeature> m_aggregate;         //! list of aggregated features

};

template<class Array>
class CInitFrameAggregator: public CDescriptorAggregator<Array> {

public:

    //! Constructor.
    CInitFrameAggregator(CTracker* tracker, const char* name):CDescriptorAggregator<Array>::CDescriptorAggregator(tracker,name){};

private:

    //! Gets the descriptor from the first feature in the tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    using CDescriptorAggregator<Array>::m_name;
    using CDescriptorAggregator<Array>::m_aggregate;

};

template<class Array>
class CSubsampleAggregator:public CDescriptorAggregator<Array> {

public:

    //! Constructor.
    CSubsampleAggregator(CTracker* tracker, const char* name, size_t n);

private:

    //! Subsamples each tracklet.
    virtual void AggregateTracklet(CTracklet* tracklet);

    size_t m_n;                          //! downsampling factor

    using CDescriptorAggregator<Array>::m_name;
    using CDescriptorAggregator<Array>::m_aggregate;

};




}
#endif // DAGG_H
