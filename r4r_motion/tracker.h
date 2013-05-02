//////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////////

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


#include <opencv2/opencv.hpp>
#include <algorithm>
#include <memory>

#include "tracklet.h"
#include "params.h"
#include "intimg.h"

namespace R4R {

/*! \brief tracker interface
 *
 *
 */
class CTracker:public std::vector<std::list<std::shared_ptr<CTracklet> > >  {

public:


	/*! \brief iterator class for CTracker objects
	 *
	 *
	 *
	 */
	class iterator {

	public:

		//! Constructor.
		iterator(CTracker* tracker);

		//! Increments iterator.
		void operator++();

		//! Increments by more than 1.
		void Advance(size_t step);

		//! Checks termination condition.
		bool operator()() { return m_t<=m_tracker->m_global_t; };

		//! Dereferencing operator.
        std::map<shared_ptr<CTracklet>,std::list<CFeature>::iterator > operator*() { return m_data; };

		//! Access to the internal counter.
		size_t GetTime() const { return m_t; }

	private:

		CTracker* m_tracker;																		//!< tracker to iterate through
        std::map<shared_ptr<CTracklet>,std::list<CFeature>::iterator > m_data;	//!< data container
		size_t m_t;																					//!< global time variable

		//! Advances the iterators that are currently alive.
		void UpdateTracklets();

		//! Looks for active tracklets at #m_t and adds them to the iterator.
		void AddTracklets();

		//! Removes dead tracks from iterator.
		void Clean();

	};

	friend class CTracker::iterator;


	//! Standard constructor.
	CTracker();

	//! Standard constructor.
    CTracker(CParameters* params);

	//! Standard destructor.
	virtual ~CTracker();

	//! Deletes all tracklets which have died from memory.
	void DeleteInvalidTracks();

	//! Returns pointer to all tracklets at given scale.
	std::list<std::shared_ptr<CTracklet> > Get(size_t s) { return at(s); };

	//! Computes the number of tracklets in the container.
	size_t Capacity();

	//! Computes the number of tracklets in the container that are still alive.
	size_t ActiveCapacity();

	/*!
	 * \brief  Initializes the tracker.
	 *
	 * \param[in] img initial image
	 *
	 */
	virtual bool Init(std::vector<cv::Mat>& pyramid) { return 0; };

	/*! \brief Executes one step of differential motion estimation.
	 *
	 * \param[in] img0 frame at t
	 * \param[in] img1 frame at t+1
	 *
	 */
	virtual bool Update(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1) { m_global_t++ ; return 0; };

	//! Draws the current features into an image.
	virtual void Draw(cv::Mat& img);

	//! Adds new features to the tracker.
	virtual bool AddTracklets(std::vector<cv::Mat>& pyramid) { return 0; };

	/*! \brief Updates all descriptors if any.
	 *
	 * \param[in] img0 frame at t
	 * \param[in] img1 frame at t+1
	 *
	 */
	virtual bool UpdateDescriptors(std::vector<cv::Mat>& pyramid) { return 0; };

	/*! \brief Marks tracks as invalid.
	 *
	 */
    virtual void Clean(std::vector<cv::Mat>& pyramid0, std::vector<cv::Mat>& pyramid1);

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
	virtual bool SaveToFile(const char* dir, const char* prefix);

	//! Searches for a tracklet with given initial feature and initial time.
    std::shared_ptr<CTracklet> SearchTracklet(CFeature x0, size_t t0);

	//! Searches for the tracklet with maximal life time.
	std::shared_ptr<CTracklet> SearchFittestTracklet();

	//! Computes the adjacency list of the covisibility graph.
	std::map<std::shared_ptr<CTracklet>,std::list<std::shared_ptr<CTracklet> > > ComputeCovisibilityGraph();

	//! Deletes descriptors attached to any of tracked features.
	void DeleteDescriptors();

	//! Sets all tracklets to active.
	void SetAllTrackletsActive();

	/*! \brief Computes integral image over locations of active features.
	 *
	 *	\details The location of active features are treated as Dirac delta functions. This routine delivers an approximation
	 *	of the indefinite integral over these functions which, upon evaluation, can be used e.g. to count the number
	 *	of active tracklets in a rectangular subdomain of the current frame. Caveat: This is just the density, meaning that
	 *	CIntegralImage::Compute() must be called separately.
	 *
	 * \returns Integral image \f$\mathcal{I}\f$. At evaluation, no swapping of indices is required, i.e., if \f$(i,j)\f$
	 * denotes a pixel where \f$i\f$ is the horizontal image location and \f$j\f$ the vertical one, the integral bound is
	 * passed in that same order \f$\mathcal{I}(i,j)\f$.
	 *
	 */
	CIntegralImage<size_t> ComputeFeatureDensity(size_t width, size_t height, size_t s);

	//! Returns the set of parameters.
    CParameters GetParameters() { return *m_params; };

	//! Draws active tracklets into an image.
	void DrawTails(cv::Mat& img, size_t length);

    //! Returns the number of tracks still alive.
    size_t GetNumberOfActiveTracks() { return m_n_active_tracks; }

	//! Lets the user select an initial bounding box.
	static cv::Rect GetManualBoundingBox(cv::Mat& img);

	//! Callback routine for selecting a bounding box.
	static void OnMouseSelectBoundingBox(int event, int x, int y, int flags, void* params);

    //! Access to the global time.
    size_t GetTime() { return m_global_t; };

    /*!
     * \brief Adds a new tracklet.
     *
     * \details Creates a new tracklet on the heap and adds a pointer to the tracker and returns
     * this pointer. Opposed to CTracklet, here, dynamic memory allocation happens within the
     * class. The destructor CTracker::~CTracker() de-allocates all tracklets that have been
     * associated with the tracker.
     *
     * \param[in] x pointer to a dynamically allocated CFeatureDescriptor object
     *
     */
    std::shared_ptr<CTracklet> AddTracklet(CFeature x);

protected:

    CParameters* m_params;					//!< container for user-defined parameters
	size_t m_global_t;						//!< global time variable
    size_t m_n_active_tracks;               //!< number of active tracks



	//! Adds a new tracklet to the pool.
	void AddTracklet(std::shared_ptr<CTracklet> tracklet);

};

}

#endif /* TRACKER_H_ */


