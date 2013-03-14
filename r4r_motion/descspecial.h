#ifndef DESCSPECIAL_H
#define DESCSPECIAL_H

#include "descriptor.h"
#include "rect.h"

#ifdef HAVE_FFTW
#include <fftw3.h>

namespace R4R {

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

    //! Paraemtrized constructor.
    CHistogramOfGradients(CRectangle<double> roi, size_t no_cells, size_t cell_size, size_t no_bins, bool osigned, double sigma);

    //! \copydoc CDescriptor::Compute(cv::Mat&)
    bool Compute(cv::Mat& img);

    //! Normalize.
    void Normalize(size_t method, double alpha);

private:

    size_t m_no_cells;              //! number of cells (per dimension)
    size_t m_cell_size;             //! number of pixels in each cell (per dimension)
    size_t m_no_bins;               //! number of histogram bins
    bool m_osigned;                  //! flag indicating whether orientations are treated as signed or unsigned
    double m_sigma;                 //! variance of Gaussian windowing function, if zero no windowing is performed

};


}


#endif // DESCSPECIAL_H
