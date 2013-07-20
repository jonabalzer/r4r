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

#ifndef R4RIMAGE_H
#define R4RIMAGE_H

#ifdef QT_GUI_LIB
#include <QImage>
#endif

#include "types.h"
#include "rect.h"

namespace R4R {

/*! \brief R4R's own gray value image class
 *
 *
 *
 */
class CImage: public CDenseArray<unsigned char> {

public:

    //! Constructor.
    CImage():CDenseArray<unsigned char>(){}

    //! Constructor.
    CImage(size_t w, size_t h):CDenseArray<unsigned char>(h,w){}

};


/*! \brief R4R's own gray value image class
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

    //! Access using bilinear interpolation.
    template<typename T> CVector<T,3> Get(const CVector<T,2>& p);

    //! Compute gradient with centered differences.
    template<typename T> CVector<T,3> Gradient(const CVector<T,2>& p, bool dir);

    //! Get the width of the image.
    size_t Width() { return NCols(); }

    //! Get the height of the image.
    size_t Height() { return NRows(); }

    //! Creates a deep copy of a region of interest.
    //CRGBImage Clone(CRectangle<int> roi);

    //! Creates a deep copy of a region of interest using bilinear interpolation.
    //template<typename T> CRGBImage Clone(CRectangle<T> roi);




private:

};

}

#endif // IMAGE_H
