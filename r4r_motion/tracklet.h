/*
 * tracklet.h
 *
 *  Created on: Mar 16, 2012
 *      Author: jbalzer
 */

#ifndef TRACKLET_H_
#define TRACKLET_H_

#include <list>
#include <set>

#include "feature.h"
#include <iostream>
#include <memory>

namespace R4R {

/*! \brief feature trajectory
 *
 *
 *
 */
class CTracklet:public std::list<CFeature> {

	friend class CTracker;

public:

	//! Standard constructor.
	CTracklet();

	//! Constructor.
    CTracklet(size_t t0, size_t s, CFeature x0);

	/*! \brief Deletes all features in the time series from memory.
	 *
	 * \details To ensure polymorphic behavior, we can only store pointers to feature objects
	 * in the tracklet container. Features must be dynamically allocated from outside
	 * of #CTracklet, and corresponding pointers passed. This function deletes them to
	 * avoid memory leaks.
	 *
	 */
	virtual ~CTracklet();

	//! Tests whether two tracklets are equal by looking at their initial position.
	bool operator!=(CTracklet& tracklet) { return GetHash()!=tracklet.GetHash(); };

	//! Directly updates the state without any filtering.
    void Update(CFeature x);

    //! Provides access to the current state.
    CFeature& GetLatestState() { return back(); };

	//! Provides access to the current feature position.
    vec GetLatestLocation() { return back().GetLocation(); };

	//! Provides access to previous feature position.
	vec GetPastLocation(size_t steps);

	//! Provides access to previous feature position w.r.t. to the native scale.
	vec GetPastLocationAtNativeScale(size_t steps);

	//! Provides access to the current feature position w.r.t. to the native scale.
    vec GetLatestLocationAtNativeScale() { return back().GetLocationAtNativeScale(); };

	//! Extrapolates the current translational speed from the feature.
	vec GetLatestVelocity();

	//! Draws the trajectory into a single image.
	void Draw(cv::Mat& img, size_t length);

	//! Streams the tracklet to a given output.
	friend std::ostream& operator<<(std::ostream& os, CTracklet& x);

	//! Writes the tracklet to a file.
	bool SaveToFile(const char* filename);

	/*! \brief Generates a hash key for the tracklet.
	 *
	 * \details The hash is built from two quantities that should identify each tracklet uniquely:
	 * - creation time
	 * - initial x and y position
	 */
	std::string GetHash();

	//! Deletes descriptors attached to any of the features in the tracklet.
	void DeleteDescriptors();

	//! Gets the status.
	bool GetStatus() { return m_status; };

	//! Sets the status.
	void SetStatus(bool status) { m_status = status; };

	/*! Gets the creation time.
	 *
	 *\details Note that there is no member function to change #m_to as it is uniquely defined during reconstruction.
	 *
	 */
	size_t GetCreationTime() { return m_t0; };

	//! Returns (initial) scale of the tracklet.
	size_t GetScale() { return m_scale; };

	//! Computes the standard deviation of feature locations over time.
	vec ComputeVariance();

protected:

	size_t m_t0;														//!< creation time
	size_t m_scale;														//!< scale
	bool m_status;														//!< status


};


/*! \brief tracklets for tracking on the selection tree
 *
 *
 *
 */
class CSTTracklet:public CTracklet {

	friend class CTST;

public:

	//! Standard constructor.
	CSTTracklet():CTracklet(),m_children(), m_orphan(true){};

	//! Constructor.
    CSTTracklet(size_t t0, size_t s, CFeature x0):CTracklet(t0,s,x0),m_children(),m_orphan(true){};

protected:

	std::set<shared_ptr<CSTTracklet> > m_children;			//!< pointer to children tracklets
	bool m_orphan;											//!< true if tracklet is not part of the tree

};

}


#endif /* TRACKLET_H_ */

