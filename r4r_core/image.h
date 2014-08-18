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

#ifndef R4RIMAGE_H
#define R4RIMAGE_H

#ifdef QT_GUI_LIB
#include <QImage>
#endif

#include <iostream>

#include "darray.h"
#include "rect.h"

namespace R4R {

/*! \brief R4R's own gray value image class
 *
 *
 *
 */
class CGrayValueImage: public CDenseArray<unsigned char> {

public:

    //! Constructor.
    CGrayValueImage():CDenseArray<unsigned char>(){}

    //! Constructor.
    CGrayValueImage(size_t w, size_t h):CDenseArray<unsigned char>(h,w){}

    //! Constructor.
    template<class Array>
    CGrayValueImage(const Array& img);

#ifdef QT_GUI_LIB
    //! Construct from QT image.
    //CGrayValueImage(const QImage& img);

    //! Cast to a QT image.
    operator QImage() const;
#endif

    //! Get the width of the image.
    size_t Width() { return NCols(); }

    //! Get the height of the image.
    size_t Height() { return NRows(); }

    //! Compute gradient with centered differences.
    std::vector<double> Gradient(const vec2& p) { return CDenseArray::Gradient<double>(p); }

    //! Compute gradient with centered differences.
    CVector<short,2> Gradient(size_t x, size_t y);

};


/*! \brief R4R's own RGB image class
 *
 *
 *
 */
class CRGBImage: public CDenseArray<rgb> {

public:

    //! Constructor.
    CRGBImage();

    //! \copydoc CDenseArray(size_t,size_t,std::shared_ptr<T>)
    CRGBImage(size_t w, size_t h, std::shared_ptr<rgb> data):CDenseArray<rgb>(h,w,data) {}

    //! Constructor.
    CRGBImage(size_t w, size_t h):CDenseArray<rgb>::CDenseArray(h,w){}

    //! Constructor.
    CRGBImage(size_t w, size_t h, unsigned char* data);

#ifdef QT_GUI_LIB
    //! Construct from QT image.
    CRGBImage(const QImage& img);

    //! Cast to a QT image.
    operator QImage() const;
#endif

    //! Get the width of the image.
    size_t Width() { return NCols(); }

    //! Get the height of the image.
    size_t Height() { return NRows(); }

    //! Compute gradient with centered differences.
    std::vector<vec3> Gradient(const vec2& p) const { return CDenseArray::Gradient<vec3>(p); }

    //! Compute gradient with centered differences.
    std::vector<vec3> Gradient(size_t i, size_t j) const { return CDenseArray::Gradient<vec3>(i,j); }

    //! Compute gradient with centered differences.
    //CVector<short,2> Gradient(size_t x, size_t y);

    //! Creates a deep copy of a region of interest.
    //CRGBImage Clone(CRectangle<int> roi);

    //! Creates a deep copy of a region of interest using bilinear interpolation.
    //template<typename T> CRGBImage Clone(CRectangle<T> roi);

private:

};

/*! \brief integral image
 *
 *
 *
 */
template<typename T>
class CIntImage: public CDenseArray<T> {

public:

    //! Constructor.
    CIntImage(size_t width, size_t height);

    //! Computes the integral.
    void Compute();

    //! Image width.
    size_t Width() { return this->NCols(); }

    //! Image height.
    size_t Height() { return this->NRows(); }

    /*! \brief Increases the density at a point.
     *
     * \param[in] x point in the image plane
     * \param[in] val amount of mass to add
     *
     * \details If the point falls between grid cells, the mass is distributed between the neighboring
     * vertices according to the individual areas of the cell partition.
     */
    template<typename U> void AddMass(const CVector<U,2>& x, T val) {

        U dx, dy;
        dx = x.Get(0) - floor(x.Get(0));
        dy = x.Get(1) - floor(x.Get(1));

        size_t i, j;
        i = size_t(x.Get(1));
        j = size_t(x.Get(0));

        this->operator ()(i,j) += T(dx*dy*val);
        this->operator ()(i+1,j) += T(dx*(1-dy)*val);
        this->operator ()(i,j+1) += T(dy*(1-dx)*val);
        this->operator ()(i+1,j+1) += T((1-dx)*(1-dy)*val);

    }

    /*! \brief Increases the density at a point.
     *
     * \param[in] x point in the image plane
     * \param[in] val amount of mass to add
     *
     */
    void AddMass(size_t i, size_t j, T val);

    /*! \brief Evaluates the integral image at corners of a rectangular window around a location.
     *
     * \param[in] x center of the integration domain
     * \param[in] hsize half-size of the integration domain
     *
     * \details \f$x\f$ is the center pixel. Non-integral locations will be interpolated bi-linearly.
     *
     */
    template<typename U,typename V> U Evaluate(const CVector<V,2>& x, const CVector<V,2>& hsize) const {

        CVector<V,2> tl = this->template ProjectToBoundary<V>(x - hsize);
        CVector<V,2> br = this->template ProjectToBoundary<V>(x + hsize);
        CVector<V,2> bl = { tl.Get(0), br.Get(1) };
        CVector<V,2> tr = { br.Get(0), tl.Get(1) };

        return this->template Get<U>(tl) + this->template Get<U>(br) - this->template Get<U>(bl) - this->template Get<U>(tr);

    }

    /*! \brief Evaluates the integral image at corners of a rectangular window around a location.
     *
     * \param[in] x center of the integration domain
     * \param[in] hsize half-size of the integration domain
     *
     * \details \f$x\f$ is the center pixel. Non-integral locations will be rounded down.
     *
     */
    template<typename U> T EvaluateApproximately(const CVector<U,2>& x, const CVector<U,2>& hsize) const {

        CVector<U,2> tl = this->ProjectToBoundary(x - hsize);
        CVector<U,2> br = this->ProjectToBoundary(x + hsize);
        CVector<U,2> bl = { tl.Get(0), br.Get(1) };
        CVector<U,2> tr = { br.Get(0), tl.Get(1) };

        return this->Get(tl) + this->Get(br) - this->Get(bl) - this->Get(tr);

    }

private:


};





}

#endif // IMAGE_H
