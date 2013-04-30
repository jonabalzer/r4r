#ifndef DESCSPECIAL_H
#define DESCSPECIAL_H

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
    bool Compute(cv::Mat& img);

};

#endif // HAVE_FFTW

class CHistogramOfGradients:public CNeighborhoodDescriptor<CRectangle<double>,vecf> {

public:

    //! Constructor.
    CHistogramOfGradients(CRectangle<double> roi);

    //! Parametrized constructor.
    CHistogramOfGradients(CRectangle<double> roi, size_t no_cells, size_t cell_size, size_t no_bins, bool osigned, double sigma);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(cv::Mat& img);

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
    CFMHoGDescriptor(CRectangle<double> roi, size_t hsize, size_t nbins);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(cv::Mat& img);

private:

    size_t m_hsize;                 //! size of the neighboorhood in which the descriptor is computed
    size_t m_no_bins;               //! number of histogram bins

};

#endif // HAVE_FFTW



}


#endif // DESCSPECIAL_H
