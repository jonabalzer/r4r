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

#ifndef R4RIDDESC_H_
#define R4RIDDESC_H_

#include "descriptor.h"
#include "rect.h"

namespace R4R {

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
	bool Compute(cv::Mat& img);

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
	bool Compute(cv::Mat& img);

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
	bool Compute(cv::Mat& img);

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
    bool Compute(cv::Mat& img);

    //! Generates sampling points.
    static void GenerateSamplePoints();

    //! Distance between two BRIEF descriptors.
    double Distance(const CBRIEF& desc) const;

protected:

    static double m_x[2*BITSET_LENGTH];             //! first set of sample points
    static double m_y[2*BITSET_LENGTH];             //! second set of sample points

};


}



#endif /* IDDESC_H_ */
