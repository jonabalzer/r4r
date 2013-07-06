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

#ifndef R4RTVAL_H_
#define R4RTVAL_H_


#include <opencv2/opencv.hpp>
#include "tracker.h"


namespace R4R {

/*! \brief tracker validation and visualization
 *
 *
 *
 */
class CTrackerValidation {

public:

	//! Constructor.
	CTrackerValidation(const char* file, cv::Size size);

	//! Constructor.
	CTrackerValidation(const char* file);

	//! Constructor.
	CTrackerValidation(const char* file, double fps);

	//! Desctructor.
	~CTrackerValidation() { m_input.release(); };

	/*! \brief Writes the video #m_input to file.
	 *
	 * \param[out] file name of output file
	 * \param[in] T maximal length of output video
	 *
	 */
	bool WriteVideo(const char* file, size_t T);

	/*! Writes a video overlayed with the tracks of different trackers.
	 *
	 * \param[in] file name of output file
	 * \param[in] trackers pointers to trackers to overlay with video stram
	 * \param[in] tmax maximum number of frames to write
	 *
	 */
	bool WriteTracks(const char* file, std::vector<CTracker*> trackers, size_t tmax);

	/*! \brief Writes a RoI to a sequence of images w.r.t. to its covariant frame.
	 *
	 * \param[out] dir output directory
	 * \param[out] prefix prefix to identify the kind of tracker used
	 * \param[in] tracklet pointer to the tracklet to save
	 * \param[in] hsize half of the size of the RoI
	 *
	 */
	bool WriteTrack(const char* dir, const char* prefix, shared_ptr<CTracklet> tracklet, size_t hsize);

	/*! \brief Shows some statistics of a tracker over its life time.
	 *
	 * \details Currently being computed are:
	 * - the number of tracks that have been created in total
	 * - maximum number of frames that some tracklet survived
	 * - average number of frames a tracklet survived
	 *
	 */
	void ShowStatistics(std::vector<CTracker*> trackers);

	//! Saves the covisibility graph to file.
	bool SaveToFile(const char* filename, std::map<shared_ptr<CTracklet>,std::list<shared_ptr<CTracklet> > > graph);

private:

	cv::VideoCapture m_input;				//!< input stream
	cv::Size m_size;						//!< size of output video
	double m_fps;							//!< frame rate of output


};

}

#endif /* TVAL_H_ */
