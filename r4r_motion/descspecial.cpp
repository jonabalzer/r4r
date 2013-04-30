#include "descspecial.h"
#include "interp.h"

using namespace cv;

namespace R4R {

#ifdef HAVE_FFTW

CFourierModulusDescriptor::CFourierModulusDescriptor(CRectangle<double> roi, size_t hsize):
    CNeighborhoodDescriptor(mat(2*hsize+1,2*hsize+1),roi) {}


bool CFourierModulusDescriptor::Compute(Mat& img) {

    fftw_complex* fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*m_container.NElems());

    double dx = 2.0/(double)(m_container.NCols()-1);
    double dy = 2.0/(double)(m_container.NRows()-1);

    for(size_t i=0; i<m_container.NRows(); i++) {

        for(size_t j=0; j<m_container.NCols(); j++) {

            vec pt = m_roi.TransformFrom(-1.0+j*dx,-1.0+i*dy);

            // row major
            fft[i*m_container.NCols() + j][0] = CImageInterpolation::Bilinear(img,pt(0),pt(1));
            fft[i*m_container.NCols() + j][1] = 0;

        }

    }

    // create plan
    fftw_plan p = fftw_plan_dft_2d(m_container.NRows(),
                                   m_container.NCols(),
                                   fft,
                                   fft,
                                   FFTW_FORWARD,
                                   FFTW_ESTIMATE);

    // run it
    fftw_execute(p);

    // compute the modulus
    for(size_t i=0; i<m_container.NRows(); i++) {

        for(size_t j=0; j<m_container.NCols(); j++) {

            m_container(i,j) = sqrt(fft[i*m_container.NCols() + j][0]*fft[i*m_container.NCols() + j][0]+fft[i*m_container.NCols() + j][1]*fft[i*m_container.NCols() + j][1]);

        }

    }

    // tidy up
    fftw_destroy_plan(p);
    fftw_free(fft);

    return 0;

}

#endif

CHistogramOfGradients::CHistogramOfGradients(CRectangle<double> roi):
    CNeighborhoodDescriptor(vecf(81),roi),
    m_no_cells(3),
    m_cell_size(6),
    m_no_bins(9),
    m_osigned(false),
    m_sigma(0) {

}

CHistogramOfGradients::CHistogramOfGradients(CRectangle<double> roi, size_t no_cells, size_t cell_size, size_t no_bins, bool osigned, double sigma):
    CNeighborhoodDescriptor(vecf(no_cells*no_cells*no_bins),roi),
    m_no_cells(no_cells),
    m_cell_size(cell_size),
    m_no_bins(no_bins),
    m_osigned(osigned),
    m_sigma(sigma) {

}

bool CHistogramOfGradients::Compute(cv::Mat& img) {

    // number of pixels in patch
    size_t n = m_no_cells*m_cell_size;
    double h = 2.0/(double)(n-1);

    // size of histogram bins
    double ho = 0;
    if(m_osigned)
        ho = 2*M_PI/m_no_bins;
    else
        ho = M_PI/m_no_bins;

    // compute actual descriptor
    for(size_t i=0; i<n; i++) {

        for(size_t j=0; j<n; j++) {

            // interpolate point
            vec pt = m_roi.TransformFrom(-1.0+j*h,-1.0+i*h);

            // compute gradient
            double Ix, Iy;
            Ix = CImageInterpolation::Gradient(img,pt(0),pt(1),0);
            Iy = CImageInterpolation::Gradient(img,pt(0),pt(1),1);

            // compute signed/unsigned orientation
            double ograd = 0;

            if(m_osigned)
                ograd = atan2(Iy,Ix) + 2*M_PI;
            else
                ograd = fabs(atan2(Iy,Ix));

            // compute the bin
            size_t bin = (size_t)ograd/ho;

            // get the cell index
            size_t I, J;
            I = (size_t)(i/m_cell_size);
            J = (size_t)(j/m_cell_size);
            size_t cindex = (J*m_no_cells + I)*(m_no_bins);   // row-major

            // increment bin weighted with gradient norm
            m_container(cindex+bin) += sqrt(Ix*Ix+Iy*Iy);

        }

    }

    return 0;

}

void CHistogramOfGradients::Normalize(size_t method, double alpha) {

    float weight;

    switch(method) {

    case 1: // regularized L1-norm, BAD
        weight = m_container.Norm1() + alpha;
        break;

    case 2: // regularized L2-norm
        weight = m_container.Norm2() + alpha;
        break;

    default: // no normalization
        weight = 1;
        break;
    }

    m_container.Scale(1.0/weight);

}


#ifdef HAVE_FFTW

CFMHoGDescriptor::CFMHoGDescriptor(CRectangle<double> roi, size_t hsize, size_t nbins):
    CNeighborhoodDescriptor(vecf(nbins),roi),
    m_hsize(hsize),
    m_no_bins(nbins) {}


bool CFMHoGDescriptor::Compute(cv::Mat& img) {

    // allocate space for histogram
    fftw_complex* fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*m_container.NElems());

    // clear bins
    for(size_t i=0; i<m_container.NElems(); i++) {

        fft[i][0] = 0;
        fft[i][1] = 0;

    }

    // full patch sidelength
    size_t n = 2*m_hsize+1;
    double h = 2.0/(double)(n-1);

    // size of each histogram bin
    double ho = 2*M_PI/m_no_bins;

    // compute actual descriptor
    for(size_t i=0; i<n; i++) {

        for(size_t j=0; j<n; j++) {

            // interpolate point
            vec pt = m_roi.TransformFrom(-1.0+j*h,-1.0+i*h);

            // compute gradient
            double Ix, Iy;
            Ix = CImageInterpolation::Gradient(img,pt(0),pt(1),0);
            Iy = CImageInterpolation::Gradient(img,pt(0),pt(1),1);

            // compute signed orientation
            double ograd = atan2(Iy,Ix) + 2*M_PI;

            // compute the bin
            size_t bin = (size_t)ograd/ho;

            // increment bin weighted with gradient norm
            fft[bin][0]  += sqrt(Ix*Ix+Iy*Iy);

        }

    }

    // compute fft
    fftw_plan p = fftw_plan_dft_1d(m_container.NElems(),
                                   fft,
                                   fft,
                                   FFTW_FORWARD,
                                   FFTW_ESTIMATE);

    fftw_execute(p);

    // compute modulus
    for(size_t i=0; i<m_container.NElems(); i++)
        m_container(i) = (float)sqrt(fft[i][0]*fft[i][0]+fft[i][1]*fft[i][1]);

    // tidy up
    fftw_destroy_plan(p);
    fftw_free(fft);

    return 0;

}

#endif

}
