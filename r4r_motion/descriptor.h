/*
 * descriptor.h
 *
 *  Created on: Apr 2, 2012
 *      Author: jbalzer
 */



#ifndef DESCRIPTOR_H_
#define DESCRIPTOR_H_

#define BITSET_LENGTH 256

#include "types.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <bitset>

namespace R4R {

/*! \brief descriptor interface
 *
 *
 *
 */
class CAbstractDescriptor {

	friend class CFeature;

public:

	//! Triggers computation of the descriptor.
	virtual bool Compute(cv::Mat& img) { return 0; };

	//! Writes descriptor to a stream.
    virtual void Write(std::ofstream& os) {};

    //! Reads a descriptor from a file stream.
    virtual void Read(std::ifstream& os) {};

    //! Visualizes the descriptor.
    virtual void Draw(cv::Mat& img, cv::Scalar color) const {};  // FIXME: get rid of this

protected:

};


template <class Array>
class CDescriptor: public CAbstractDescriptor {

public:

    //! Constructor.
    CDescriptor();

    //! Constructor.
    CDescriptor(const Array& data);

    //! \copydoc CAbstractDescriptor::Compute(cv::Mat&)
    virtual bool Compute(cv::Mat& img) { return 0; };

    //! \copydoc CAbstractDescriptor::Write(std::ofstream&)
    virtual void Write(std::ofstream& ofs) { ofs << m_container; };

    //! \copydoc CAbstractDescriptor::Read(std::ifstream&)
    virtual void Read(std::ifstream& ifs) { ifs >> m_container; };

     //! Accesses the descriptor container.
    Array& Get() { return m_container; };

    //! Sets the descriptor data from the outside.
    void Set(const Array& container) { m_container = container; };

protected:

    Array m_container;

};


template <class Rect, class Array>
class CNeighborhoodDescriptor:public CDescriptor<Array> {

public:

	//! Constructor.
    CNeighborhoodDescriptor(const Rect& roi);

    //! Constructor.
    CNeighborhoodDescriptor(const Array& container, const Rect& roi);

	//! \copydoc CDescriptor::Draw(cv::Mat& img,cv::Scalar)
    virtual void Draw(cv::Mat& img, cv::Scalar color) { m_roi.Draw(img,color); };

    //! Access to the region of interest.
    Rect GetRoI() const { return m_roi; };

protected:

    Rect m_roi;										//!< region of interest

    using CDescriptor<Array>::m_container;

};

}

#endif /* DESCRIPTOR_H_ */
