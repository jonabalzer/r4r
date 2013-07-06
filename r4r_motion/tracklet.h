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

#ifndef R4RTRACKLET_H_
#define R4RTRACKLET_H_

#include <list>
#include <set>

#ifdef QT_GUI_LIB
#include <QImage>
#include <QPainter>
#include <QPen>
#endif

#include "feature.h"
#include "types.h"
#include <iostream>
#include <memory>

namespace R4R {

typedef CInterestPoint<float,2> imfeature;

/*! \brief feature trajectory
 *
 *
 *
 */
class CTracklet:public std::list<imfeature> {

	friend class CTracker;

public:

	//! Standard constructor.
	CTracklet();

	//! Constructor.
    CTracklet(size_t t0, imfeature x0);

	//! Tests whether two tracklets are equal by looking at their initial position.
    bool operator!=(CTracklet& tracklet) { return GetHash()!=tracklet.GetHash(); }

	//! Directly updates the state without any filtering.
    void Update(imfeature x);

    //! Provides access to the current state.
    imfeature& GetLatestState() { return back(); }

	//! Provides access to the current feature position.
    vec2f GetLatestLocation() { return back().GetLocation(); }

	//! Provides access to previous feature position.
    vec2f GetPastLocation(size_t steps);

	//! Provides access to previous feature position w.r.t. to the native scale.
    vec2f GetPastLocationAtNativeScale(size_t steps);

	//! Provides access to the current feature position w.r.t. to the native scale.
    vec2f GetLatestLocationAtNativeScale() { return back().GetLocationAtNativeScale(); }

	//! Extrapolates the current translational speed from the feature.
    vec2f GetLatestVelocity();

    //! Returns the current scale.
    float GetScale() { return back().GetScale(); }

	//! Streams the tracklet to a given output.
	friend std::ostream& operator<<(std::ostream& os, CTracklet& x);

	/*! \brief Generates a hash key for the tracklet.
	 *
	 * \details The hash is built from two quantities that should identify each tracklet uniquely:
	 * - creation time
	 * - initial x and y position
	 */
	std::string GetHash();

	//! Gets the status.
    bool GetStatus() { return m_status; }

    //! Sets the status flag.
    void SetStatus(bool status) { m_status = status; }

	/*! Gets the creation time.
	 *
     *\details Note that there is no member function to change #m_t as it is uniquely defined during reconstruction.
	 *
	 */
    size_t GetCreationTime() { return m_t0; }

#ifdef QT_GUI_LIB

    static const Qt::GlobalColor COLORS[10];

    //! Draws the trajectory into a single image.
    void Draw(QImage& img, size_t length);

#endif

protected:

    size_t m_t0;            //!< creation time
    bool m_status;			//!< status

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
    CSTTracklet():CTracklet(),m_children(), m_orphan(true){}

	//! Constructor.
    CSTTracklet(size_t t0, imfeature x0):CTracklet(t0,x0),m_children(),m_orphan(true){}

protected:

	std::set<shared_ptr<CSTTracklet> > m_children;			//!< pointer to children tracklets
	bool m_orphan;											//!< true if tracklet is not part of the tree

};

} // end of namespace

#endif /* TRACKLET_H_ */

