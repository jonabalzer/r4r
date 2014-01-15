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
#include <iostream>
#include <memory>
#include <iterator>

#ifdef QT_GUI_LIB
#include <QImage>
#include <QPainter>
#include <QPen>
#endif

#include "feature.h"
#include "types.h"

#define MAX_TRACKLET_LENGTH 50

namespace R4R {

typedef CInterestPoint<float,2> imfeature;

/*! \brief feature trajectory
 *
 *
 */
class CTracklet {

public:

    //! Inhibits empty tracklets.
    CTracklet() = delete;

    //! Constructor.
    explicit CTracklet(size_t t0, const imfeature& x0, size_t maxlength = MAX_TRACKLET_LENGTH);

    //! Directly updates the state without any filtering.
    void Update(const imfeature& x);

    /*! \brief Provides access to previous feature position.
     *
     * \param[in] steps number of steps to go back in time
     *
     */
    const vec2f& GetPastLocation(size_t steps) const;

    //! Provides access to previous feature position w.r.t. to the native scale.
    vec2f GetPastLocationAtNativeScale(size_t steps) const;

    //! Streams the tracklet to a given output.
    friend  std::ostream& operator <<(std::ostream& os, const CTracklet& x);

    //! Returns the lifetime of the tracklet.
    size_t Size() { return m_lifetime; }

    /*! \brief Tests whether two tracklets are equal by looking at their initial position.
     *
     * This can always be done as there is no way to construct an empy tracklet.
     *
     */
    bool operator!=(const CTracklet& tracklet) const { return GetHash()!=tracklet.GetHash(); }

    //! Provides access to the current state.
    imfeature& GetLatestState() { return m_data.at(m_cursor); }

    //! Provides access to the current state.
    const imfeature& GetLatestState() const { return m_data.at(m_cursor); }

    //! Provides access to the current feature position.
    const vec2f& GetLatestLocation() const { return m_data.at(m_cursor).GetLocation(); }

    //! Provides access to the current feature position w.r.t. to the native scale.
    vec2f GetLatestLocationAtNativeScale() const { return m_data.at(m_cursor).GetLocationAtNativeScale(); }

    //! Returns the current scale.
    float GetScale() const { return m_data.at(m_cursor).GetScale(); }

    //! Gets the status.
    bool GetStatus() const { return m_status; }

    //! Sets the status flag.
    void SetStatus(bool status) { m_status = status; }

    /*! \brief Gets the creation time.
     *
     * Note that there is no member function to change #m_t as it is uniquely defined during reconstruction.
     *
     */
    size_t GetCreationTime() const { return m_t0; }

    //! Returns tracklet hash code.
    std::string GetHash() const { return m_hash; }

    //! Read-access to the data.
    const vector<imfeature>& GetData() { return m_data; }

    /*! \brief Generates a hash key for the tracklet.
     *
     * The hash is built from two quantities that should identify each tracklet uniquely:
     * - creation time
     * - initial location
     */
     static std::string GenerateHash(size_t t0, const imfeature& x);

     /*! \brief Deletes unused ring buffer space.
      *
      * This method is useful for exporting entire trajectories. In this case, choose the buffer size
      * larger than necessary, run the tracker, and compress the tracklets before saving them to disk.
      * This also reverses the storage order.
      *
      */
     void CompressAndReverse();

#ifdef QT_GUI_LIB

    static const Qt::GlobalColor COLORS[10];

    //! Draws the trajectory into a single image.
    void Draw(QImage& img, size_t length) const;

#endif

private:

    vector<imfeature> m_data;       //!< vector which holds the data
    size_t m_t0;                    //!< creation time
    bool m_status;                  //!< status
    std::string m_hash;             //!< hash key for tracklet
    int m_cursor;                   //!< iterator pointing to the last feature
    size_t m_lifetime;              //!< how long is the tracklet alive

};



/*! \brief tracklets for tracking on the selection tree
 *
 * Alhtough TST is not part of R4R anymore, this class is kept because it might be useful
 * anywhere where tracklets have to refer to each other.
 *
 */
class CSTTracklet:public CTracklet {

    //friend class CTST;

public:

	//! Standard constructor.
    CSTTracklet() = delete;

	//! Constructor.
    CSTTracklet(size_t t0, imfeature x0, size_t maxlength = MAX_TRACKLET_LENGTH):CTracklet(t0,x0,maxlength),m_children(),m_orphan(true){}

protected:

	std::set<shared_ptr<CSTTracklet> > m_children;			//!< pointer to children tracklets
	bool m_orphan;											//!< true if tracklet is not part of the tree

};

} // end of namespace

#endif /* TRACKLET_H_ */

