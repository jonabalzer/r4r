/*
 * feature.h
 *
 *  Created on: Apr 2, 2012
 *      Author: jbalzer
 */

#ifndef FEATURE_H_
#define FEATURE_H_

#include <opencv2/opencv.hpp>
#include "rect.h"
#include "darray.h"
#include <stdio.h>
#include <list>
#include <iostream>
#include <fstream>
#include <memory>

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

	//! Orders two image points.
	bool operator<(CFeature& x);

	//! Compares two image points.
    bool operator!=(CFeature& x);

	//! Compares two image points.
	bool operator==(CFeature& x) {  return !(*this!=x); };

	//! Returns the scale.
	size_t GetScale() { return m_scale; }

	//! Returns the strength/quality of the feature.
	double GetQuality() { return m_quality; }

	//! Sets the strength/quality of the feature.
	void SetQuality(double quality) { m_quality = quality; };

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
	size_t NoDescriptors() { return m_descriptors.size(); };

	//! Checks whether a descriptor of a given name exists.
	bool HasDescriptor(const char* name);

	//! Returns feature location w.r.t. to inherent scale.
	vec GetLocation() { return m_location; };

    //! Returns feature location w.r.t. to native scale.
	vec GetLocationAtNativeScale();

    //! Access to the descriptor container.
    std::map<string,shared_ptr<CAbstractDescriptor> >& GetDescriptors() { return m_descriptors; };

    //! Looks for descriptor with a specified name.
	shared_ptr<CAbstractDescriptor> GetDescriptor(const char* name);

	//! Looks for a descriptor a specified position.
	shared_ptr<CAbstractDescriptor> GetDescriptor(size_t no);

    //! Looks for a descriptor a specified position.
    string GetDescriptorName(size_t no);

    //! Writes a set of features to file.
    static bool SaveToFile(const char* filename, std::list<CFeature>& features);

    /*! \brief Reads a set of features from file.
     * \param[in] filename file name
     * \param[out] featurs list of features
     * \param[in] type precision of the container
     */
    static bool OpenFromFile(const char* filename, std::list<CFeature>& features, int type = F4S);

    /*! \brief Writes a set of features to file.
     *
     * \param[in] filename file name
     * \param[in] features list of features with float descriptors attached
     * \param[in] name name of the descriptor to save
     * \param[in] comment comment
     * \param[in] type data type of descriptor
     * \param[in] t0 creation time if the descriptor came from a tracklet
     * \details This is for fast input/output of massive amounts of data. Only float containers
     * can be handled.
    */
    static bool SaveDescriptors(const char* filename, std::list<CFeature>& features, const char* name, const char* comment, int type = 0, size_t t0 = 0);

    /*! \brief Reads a set of features from file.
     * \param[in] filename file name
     * \param[out] headers list of headers
     * \param[out] data pointer to a contiguous block of memory holding the descriptor data
     *
     */
    static float* LoadDescriptors(const char* filename, std::vector<CDescriptorFileHeader>& headers);

private:

    vec m_location;                                                         //! feature location
    size_t m_scale;                                                         //!< scale
    double m_quality;                                                       //!< quality measure, e.g., for storing motion estimation residual
    std::map<string,shared_ptr<CAbstractDescriptor> > m_descriptors;        //!< attached descriptors

};



}

#endif /* FEATURE_H_ */
