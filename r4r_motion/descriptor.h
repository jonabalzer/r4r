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

#ifndef R4RDESCRIPTOR_H_
#define R4RDESCRIPTOR_H_

#define BITSET_LENGTH 256

#include "types.h"
#include "feature.h"
#include <opencv2/opencv.hpp>
#include <iostream>
#include <fstream>
#include <bitset>
#include <string>

namespace R4R {

/*! \brief descriptor interface
 *
 *
 *
 */
class CAbstractDescriptor {

public:

	//! Triggers computation of the descriptor.
    virtual bool Compute(cv::Mat& img) = 0;

	//! Writes descriptor to a stream.
    virtual void Write(std::ofstream& os) = 0;

    //! Reads a descriptor from a file stream.
    virtual void Read(std::ifstream& os) = 0;

    //! Number of rows of container.
    virtual size_t NRows() = 0;

    //! Number of columns of container.
    virtual size_t NCols() = 0;

    //! Number of elements in the container.
    virtual size_t NElems() = 0;

    //! Returns the data type.
    virtual ETYPE GetType() = 0;

    //! Return pointer to the descriptor data.
    virtual void* GetData() = 0;

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
    virtual bool Compute(cv::Mat& img) { return 0; }

    //! \copydoc CAbstractDescriptor::Write(std::ofstream&)
    virtual void Write(std::ofstream& ofs) { ofs << m_container; }

    //! \copydoc CAbstractDescriptor::Read(std::ifstream&)
    virtual void Read(std::ifstream& ifs) { ifs >> m_container; }

     //! Accesses the descriptor container.
    Array& Get() { return m_container; }

    //! Sets the descriptor data from the outside.
    void Set(const Array& container) { m_container = container; }

    //! \copydoc CAbstractDescriptor::NRows()
    virtual size_t NRows() { return m_container.NRows(); }

    //! \copydoc CAbstractDescriptor::NCols()
    virtual size_t NCols() { return m_container.NCols(); }

    //! \copydoc CAbstractDescriptor::NElems()
    virtual size_t NElems() { return m_container.NElems(); }

    //! \copydoc CAbstractDescriptor::GetType()
    virtual ETYPE GetType() { return m_container.GetType(); }

    //! \copydoc CAbstractDescriptor::GetData()
    virtual void* GetData() { return m_container.Data().get(); }


protected:

    Array m_container;

private:

    // prevent copying of data
    CDescriptor(const CDescriptor<Array>& desc);
    CDescriptor<Array> operator=(const CDescriptor<Array>& desc);

};


template <class Rect, class Array>
class CNeighborhoodDescriptor:public CDescriptor<Array> {

public:

	//! Constructor.
    CNeighborhoodDescriptor(const Rect& roi);

    //! Constructor.
    CNeighborhoodDescriptor(const Array& container, const Rect& roi);

	//! \copydoc CDescriptor::Draw(cv::Mat& img,cv::Scalar)
    virtual void Draw(cv::Mat& img, cv::Scalar color) { m_roi.Draw(img,color); }

    //! Access to the region of interest.
    Rect GetRoI() const { return m_roi; }

protected:

    Rect m_roi;										//!< region of interest

    using CDescriptor<Array>::m_container;

};

class CDescriptorFileHeader {

    friend class CFeature;

public:

    //! Constructor.
    CDescriptorFileHeader();

    //! Constructor.
    CDescriptorFileHeader(const CFeature& feature);

    //! Sets the comment.
    void SetComment(const char* comment) { m_comment = std::string(comment); }

    //! Sets the descriptor name, size, and type.
    bool SetDescriptor(CFeature& feature, const char* name);

    //! Writes header to a file stream.
    friend std::ostream& operator<<(std::ostream& os, CDescriptorFileHeader& x);

    //! Reads header from a file stream.
    friend std::istream& operator>>(std::istream& is, CDescriptorFileHeader& x);

    //! Get the number of elements in the descriptor.
    size_t NElems() { return m_size[0]*m_size[1]; }

    //! Get the number of rows.
    size_t NRows() { return m_size[0]; }

    //! Get the number of rows.
    size_t NCols() { return m_size[1]; }

    //! Access to comments.
    std::string GetComment() { return m_comment; }

    //! Access to location.
    void GetLocation(float& u, float& v) { u = m_location[0]; v = m_location[1]; }

    //! Access to scale.
    size_t GetScale() { return m_scale; }

private:

    float m_location[2];                        //! location of the feature where the descriptor was computed
    size_t m_scale;                             //! scale of the mother feature
    std::string m_name;                         //! name of the descriptor
    size_t m_size[2];                           //! size of the container
    ETYPE m_type;                               //! data type
    std::string m_comment;                      //! comments

};


}

#endif /* DESCRIPTOR_H_ */
