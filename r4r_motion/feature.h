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

#ifndef R4RFEATURE_H_
#define R4RFEATURE_H_

#include <stdio.h>
#include <list>
#include <iostream>
#include <fstream>
#include <memory>

#include <opencv2/opencv.hpp>

#include "rect.h"
#include "darray.h"


using namespace std;

namespace R4R {

class CAbstractDescriptor;
class CDescriptorFileHeader;

/*! \brief feature
 *
 * \details
 *
 * \TODO Rename this class. It actually represents interest points, irrespective of dimension. What about features in
 * n dimensions? Implement generalization with polymorphism or generically?
 *
*/
class CFeature {

    friend class CAbstractDescriptor;
    friend class CDescriptorFileHeader;

public:

	static const cv::Scalar COLORS[12];						//!< color table for drawing features at different scales

	//! Standard constructor.
	CFeature();

	//! Parametrized constructor for 2d features.
	CFeature(vec loc, size_t scale, double quality = 0);

	//! Parametrized constructor for 2d features.
	CFeature(double x, double y, size_t scale, double quality = 0);

	//! Destructor.
	~CFeature();

    //! Copy constructor.
    CFeature(const CFeature& x);

	//! Orders two image points.
    bool operator<(CFeature& x);

	//! Compares two image points.
    bool operator!=(CFeature& x);

	//! Compares two image points.
    bool operator==(CFeature& x) {  return !(*this!=x); }

	//! Returns the scale.
	size_t GetScale() { return m_scale; }

	//! Returns the strength/quality of the feature.
	double GetQuality() { return m_quality; }

	//! Sets the strength/quality of the feature.
    void SetQuality(double quality) { m_quality = quality; }

    //! Writes a feature to an output stream.
    friend std::ostream& operator<<(std::ostream& os, CFeature& x);

    //! Writes a feature to a file stream.
    friend std::ofstream& operator<<(std::ofstream& os, CFeature& x);

    //! Reads a from a file stream.
    friend std::ifstream& operator>>(std::ifstream& is, CFeature& x);

	//! Draws the feature into the corresponding image at native scale.
	void Draw(cv::Mat& img);

    //! Draws the feature and into the corresponding image native scale.
	void Draw(cv::Mat& img, cv::Scalar color, size_t size = 5);

	//! Adds a descriptor to the feature.
	void AttachDescriptor(const char* id, shared_ptr<CAbstractDescriptor> descriptor);

	//! Deletes all attached descriptors from memory.
	void DeleteDescriptors();

	//! Returns the number of descriptors attached to the feature.
    size_t NoDescriptors() { return m_descriptors.size(); }

	//! Checks whether a descriptor of a given name exists.
	bool HasDescriptor(const char* name);

	//! Returns feature location w.r.t. to inherent scale.
    vec GetLocation() { return m_location; }

    //! Returns feature location w.r.t. to native scale.
    vec GetLocationAtNativeScale();

    //! Access to the descriptor container.
    std::map<string,shared_ptr<CAbstractDescriptor> >& GetDescriptors() { return m_descriptors; }

    //! Looks for descriptor with a specified name.
	shared_ptr<CAbstractDescriptor> GetDescriptor(const char* name);

    //! Looks for a descriptor at specified position.
	shared_ptr<CAbstractDescriptor> GetDescriptor(size_t no);

    //! Looks for name of descriptor at specified position.
    string GetDescriptorName(size_t no);

    //! Writes a set of features to file.
    static bool SaveToFile(const char* filename, std::list<CFeature>& features);

    /*! \brief Reads a set of features from file.
     * \param[in] filename file name
     * \param[out] featurs list of features
     * \param[in] type precision of the container
     */
    static bool OpenFromFile(const char* filename, std::list<CFeature>& features);

    /*! \brief Reads a set of features from file into a contiguous memory block.
     * \param[in] filename file name
     * \param[out] headers list of headers
     * \param[out] data pointer to a contiguous block of memory holding the descriptor data
     * \param[out] type data type to cast the void pointer to
     *
     */
    static void* LoadDescriptors(const char* filename, std::vector<CDescriptorFileHeader>& headers, ETYPE& type, bool preview = false);

//    /*! \brief Writes a set of features to file.
//     *
//     * \param[in] filename file name
//     * \param[in] features list of features with float descriptors attached
//     * \param[in] name name of the descriptor to save
//     * \param[in] comment comment
//     * \param[in] type data type of descriptor
//     * \param[in] t0 creation time if the descriptor came from a tracklet
//     * \details This is for fast input/output of massive amounts of data. Only float containers
//     * can be handled.
//    */
//    static bool SaveDescriptors(const char* filename, std::list<CFeature>& features, const char* name, const char* comment, ETYPE type = F4S, size_t t0 = 0);

    /*! \brief Reads a set of features from file.
     * \param[in] filename file name
     * \param[out] headers list of headers
     * \param[out] data pointer to a contiguous block of memory holding the descriptor data
     *
     */
    //static float* LoadDescriptors(const char* filename, std::vector<CDescriptorFileHeader>& headers);

    /*! \brief Reads a set of features from file.
     * \param[in] filename file name
     * \param[out] headers list of headers
     * \param[out] data pointer to a contiguous block of memory holding the descriptor data
     * \param[out] type data type to cast the void pointer to
     *
     */
    //static bool LoadDescriptors(const char* filename, std::vector<CDescriptorFileHeader>& headers);//, void* data, ETYPE& type);

private:

    vec m_location;                                                         //! feature location
    size_t m_scale;                                                         //!< scale
    double m_quality;                                                       //!< quality measure, e.g., for storing motion estimation residual
    std::map<string,shared_ptr<CAbstractDescriptor> > m_descriptors;        //!< attached descriptors

};



}

#endif /* FEATURE_H_ */
