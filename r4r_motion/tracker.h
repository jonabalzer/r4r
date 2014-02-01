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

/*! \mainpage
 *  \author    J. Balzer
 *  \version   0.1
 *  \date      2013-05-01
 *
 *
 *
 *
 */

#ifndef R4RTRACKER_H_
#define R4RTRACKER_H_

#ifdef QT_GUI_LIB
#include <QImage>
#endif

#include <opencv2/opencv.hpp>

#include <algorithm>
#include <memory>
#include <iterator>

#include "tracklet.h"
#include "params.h"
#include "image.h"
#include "dagg.h"

namespace R4R {

/*! \brief tracker interface
 *
 *
 */

template<template<class Tracklet, class Allocator = std::allocator<Tracklet> > class TrackerContainer,template<class T, class Allocator = std::allocator<T> > class TrackletContainer>
class CTracker {

public:

	//! Standard constructor.
    CTracker();

	//! Standard constructor.
    CTracker(CParameters* params);

    //! Standard destructor.
	virtual ~CTracker();

	//! Deletes all tracklets which have died from memory.
    void DeleteInvalidTracks();

    //! Computes the number of tracklets in the container.
    size_t Capacity() { return m_data.size(); }

    /*!
     * \brief Adds a new tracklet.
     *
     * Creates a new tracklet on the heap and adds a pointer to the tracker and returns
     * this pointer. Opposed to CTracklet, here, dynamic memory allocation happens within the
     * class. The destructor CTracker::~CTracker() de-allocates all tracklets that have been
     * associated with the tracker.
     *
     * \param[in] x feature point
     *
     */
    std::shared_ptr<CTracklet<TrackletContainer> > AddTracklet(const imfeature& x);

    //! Creates a shared pointer and adds it to the tracklet pool.
    void AddTracklet(CTracklet<TrackletContainer>* trackler);

	//! Computes the number of tracklets in the container that are still alive.
    size_t ActiveCapacity() const;

	/*!
	 * \brief  Initializes the tracker.
	 *
	 * \param[in] img initial image
	 *
	 */
    virtual bool Init(std::vector<cv::Mat>& pyramid) = 0;

	/*! \brief Executes one step of differential motion estimation.
	 *
     * \param[in] img0 frame at \f$t\f$
     * \param[in] img1 frame at \f$t+1\f$
	 *
	 */
    virtual bool Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1) = 0;

	//! Adds new features to the tracker.
    virtual bool AddTracklets(std::vector<cv::Mat>& pyramid) = 0;

	/*! \brief Updates all descriptors if any.
	 *
	 * \param[in] img0 frame at t
	 * \param[in] img1 frame at t+1
	 *
	 */
    virtual bool UpdateDescriptors(std::vector<cv::Mat>& pyramid) = 0;

	/*! \brief Marks tracks as invalid.
	 *
	 */
    virtual void Clean(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1) = 0;

	/*! \brief Saves all tracklets to file.
	 *
	 * \details A file name is created for each tracklet in the tracker automatically combining the user-defined prefix
	 * with the creation time and initial location of the feature.
	 *
	 * \param[in] dir save directory
	 * \param[in] prefix prefix that identifies the tracker in the file name
     *
     * \todo Renames this into SaveToFiles(), or make lower function part of descriptor aggregation class.
	 *
	 */
    bool SaveToFile(const char* dir, const char* prefix);

	//! Searches for a tracklet with given initial feature and initial time.
    std::shared_ptr<CTracklet<TrackletContainer> > SearchTracklet(const std::string& hash) const;

	//! Searches for the tracklet with maximal life time.
    std::shared_ptr<CTracklet<TrackletContainer> > SearchFittestTracklet() const;

	//! Computes the adjacency list of the covisibility graph.
    std::map<std::shared_ptr<CTracklet<TrackletContainer> >,std::list<std::shared_ptr<CTracklet<TrackletContainer> > > > ComputeCovisibilityGraph() const;

	//! Sets all tracklets to active.
    void SetAllTrackletsActive();

	/*! \brief Computes integral image over locations of active features.
	 *
	 *	\details The location of active features are treated as Dirac delta functions. This routine delivers an approximation
	 *	of the indefinite integral over these functions which, upon evaluation, can be used e.g. to count the number
	 *	of active tracklets in a rectangular subdomain of the current frame. Caveat: This is just the density, meaning that
	 *	CIntegralImage::Compute() must be called separately.
	 *
     * \returns
     *
	 */
    std::vector<size_t> ComputeFeatureDensity(std::vector<CIntImage<size_t> >& imgs) const;

	//! Returns the set of parameters.
    const CParameters& GetParameters() { return *m_params; }

    //! Access to the global time.
    size_t GetTime() { return m_global_t; }

    //! Read-only access to data.
    const TrackerContainer<std::shared_ptr<CTracklet<TrackletContainer> > >& GetData() { return m_data; }

    //! Adjusts the size of all tracklet buffer.
    void ResizeTracklets(size_t n);

#ifdef QT_GUI_LIB
    //! Draws active tracklets into an image.
    void Draw(QImage& img, size_t length) const;
#endif

    //! Triggers aggregation of all tracklets.
    template<class Array> std::list<imfeature> Aggregate(const CDescriptorAggregator<Array,TrackletContainer>& aggregator, const string& name) const;

protected:

    TrackerContainer<std::shared_ptr<CTracklet<TrackletContainer> > >  m_data;     //!< container holding the tracklets
    CParameters* m_params;                                                         //!< container for user-defined parameters
    size_t m_global_t;                                                             //!< global time variable

};

typedef CTracker<list,CRingBuffer> CSlidingWindowTracker;
typedef CTracker<list,list> CContinuousTracker;

} // end of namespace

#endif /* TRACKER_H_ */


