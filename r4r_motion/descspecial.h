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

#ifndef R4RDESCSPECIAL_H
#define R4RDESCSPECIAL_H

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif // HAVE_FFTW

#include "descriptor.h"
#include "rect.h"


namespace R4R {

#ifdef HAVE_FFTW

/*! \brief modulus of the FFT of the patch
 *
 *
 *
 */
class CFourierModulusDescriptor:public CNeighborhoodDescriptor<CRectangle<double>,mat> {

public:

    //! Constructor.
    CFourierModulusDescriptor(CRectangle<double> roi, size_t hsize);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

};

#endif // HAVE_FFTW

class CHistogramOfGradients:public CNeighborhoodDescriptor<CRectangle<double>,vecf> {

public:

    //! Constructor.
    CHistogramOfGradients(CRectangle<double> roi);

    //! Parametrized constructor.
    CHistogramOfGradients(CRectangle<double> roi, size_t no_cells, size_t cell_size, size_t no_bins, bool osigned, double sigma);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

    //! Normalize.
    void Normalize(size_t method, double alpha);

private:

    size_t m_no_cells;              //! number of cells (per dimension)
    size_t m_cell_size;             //! number of pixels in each cell (per dimension)
    size_t m_no_bins;               //! number of histogram bins
    bool m_osigned;                 //! flag indicating whether orientations are treated as signed or unsigned
    double m_sigma;                 //! variance of Gaussian windowing function, if zero no windowing is performed

};


#ifdef HAVE_FFTW

class CFMHoGDescriptor:public CNeighborhoodDescriptor<CRectangle<double>,vecf> {

public:

    //! Constructor.
    CFMHoGDescriptor(CRectangle<double> roi, size_t hsize = 30, size_t nbins = 128);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(const cv::Mat& img);

private:

    size_t m_hsize;                 //! size of the neighboorhood in which the descriptor is computed
    size_t m_no_bins;               //! number of histogram bins

};

#endif // HAVE_FFTW



}


#endif // DESCSPECIAL_H
