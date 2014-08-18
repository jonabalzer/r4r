//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
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
    virtual bool Compute(const cv::Mat& img) = 0;

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
    virtual bool Compute(const cv::Mat& img) { return 0; }

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

    //! Access to the region of interest.
    const Rect& GetRoI() const { return m_roi; }

protected:

    Rect m_roi;										//!< region of interest

    using CDescriptor<Array>::m_container;

};


/*! \brief descriptor that maps an image (patch) to itself
 *
 *
 *
 */
class CIdentityDescriptor:public CNeighborhoodDescriptor<CRectangle<double>,matf> {

public:

    //! Constructor.
    CIdentityDescriptor(CRectangle<double> roi, size_t method, size_t hsize = 7);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

    //! Normalizes the image patch by subtracting the mean and dividing by the standard deviation.
    void Normalize();

    //! Normalizes the image patch by dividing by the mean.
    void NormalizeWeber();

    //! Subtracts the mean from the image patch.
    void Center();

    //! Computes the correlation coefficient between two normalized images.
    double Distance(CDescriptor<matf>& desc) const;

private:

    size_t m_method;								//! flag indicating which normalization method to use


};

class CIdentityGradientDescriptor:public CNeighborhoodDescriptor<CRectangle<double>,mat> {

public:

    //! Constructor.
    CIdentityGradientDescriptor(CRectangle<double> roi, double alpha, size_t method, size_t hsize = 7);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

    /*! \brief Distance between two descriptors.
     *
     * \details Forms the scalar product of the gradient in one descriptor with the dual gradient of the other. If both
     * gradient point into the same direction (or if they both have zero norm), then the contribution to the overall sum
     * is zero.
     *
     * When both gradient fields are unit norm, the distance can be rendered rotation-invariant by first computing all
     * scalar products on the support of the descriptor then making it mean-free before accumulating the absolute value.
     *
     */
    double Distance(CDescriptor<mat>& desc) const;

protected:

    //! Normalizes all vectors which are greater than a threshold in norm.
    void Normalize();

    //! Weighting function depending on norm of gradient.
    double WeightingFunction(double gnorm);

    double m_alpha; 								//! normalization threshold
    size_t m_method;								//! flag indicating which normalization method to use


};

class CCurvatureDescriptor:public CIdentityGradientDescriptor {

public:

    //! Constructor.
    CCurvatureDescriptor(CRectangle<double> roi, double alpha, size_t method, size_t hsize = 7);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

    //! Distance between two curvature  descriptors.
    double Distance(const CCurvatureDescriptor& desc) const;

private:

    mat m_kappa;                        // container for

};


/*! \brief topological descriptor
 *
 *
 *
 */
/*class CFBDDescriptor:public CNeighborhoodDescriptor<CRectangle<double>,vec> {

public:

    //! Constructor.
    CFBDDescriptor(CRectangle<double> roi, size_t length);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(cv::Mat& img);

    //! Computes the correlation coefficient between two normalized images (not a distance in the strict sense).
    virtual double Distance(CNeighborhoodDescriptor<CRectangle<double>,vec>& desc) const;

private:

    size_t m_length;				//!< number of

};*/

/*! \brief implementation of BRIEF descriptor proposed in [Calonder2010]
 *
 *
 *
 */
class CBRIEF:public CNeighborhoodDescriptor<CRectangle<double>,CDenseVector<bool> > {

public:

    //! Constructor.
    CBRIEF(CRectangle<double> roi);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

    //! Generates sampling points.
    static void GenerateSamplePoints();

    //! Print sampling points.
    static void PrintSamplePoints();

    //! Distance between two BRIEF descriptors.
    double Distance(const CBRIEF& desc) const;

protected:

    static vec2 m_pts_0[BITSET_LENGTH];             //! first set of sample points
    static vec2 m_pts_1[BITSET_LENGTH];             //! second set of sample points

};


}

#endif /* DESCRIPTOR_H_ */
