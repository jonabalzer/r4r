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

#include "image.h"
#include <string.h>

using namespace std;

namespace R4R {

#ifdef QT_GUI_LIB
CGrayValueImage::operator QImage() const {

    QImage qimg(NCols(),NRows(),QImage::Format_RGB888);

    for(size_t i=0; i<NRows(); i++) {

        for(size_t j=0; j<NCols(); j++) {

            unsigned char val = CDenseArray<unsigned char>::Get(i,j);

            qimg.setPixel(j,i,qRgb(val,val,val));

        }

    }

    return qimg;

}
#endif

template <class Array>
CGrayValueImage::CGrayValueImage(const Array& img):
    CDenseArray<unsigned char>::CDenseArray(img.NRows(),img.NCols()) {

    for(size_t i=0; i<img.NRows(); i++) {

        for(size_t j=0; j<img.NCols(); j++)
            this->Set(i,j,(unsigned char)img.Get(i,j));

    }

}

template CGrayValueImage::CGrayValueImage<matf>(const matf& x);

#ifdef QT_GUI_LIB
template <>
CGrayValueImage::CGrayValueImage(const QImage& img):
    CDenseArray<unsigned char>::CDenseArray(img.height(),img.width()) {

    for(size_t i=0; i<NRows(); i++) {

        for(size_t j=0; j<NCols(); j++) {

            QRgb val = img.pixel(j,i);

            // find mean
            float mean = (float(qRed(val)) + float(qGreen(val)) + float(qBlue(val)))/3.0;

            // roound
            this->operator ()(i,j) = (unsigned char)(mean+0.5);

        }

    }

}
#endif

CVector<short,2> CGrayValueImage::Gradient(size_t x, size_t y) {

    if(x<1 && y<1 && x>=Width()-1 && y>=Height()-1)
        return CVector<short,2>();

        short I0x, I0y, I1x, I1y;
        I0x = short(this->Get(y,x-1));
        I0y = short(this->Get(y-1,x));
        I1x = short(this->Get(y,x+1));
        I1y = short(this->Get(y+1,x));

        short dIx, dIy;
        dIx = I1x - I0x;
        dIy = I1y - I0y;

        if(!(dIx&1)) // even number
            dIx /= 2;
        else
            dIx = (dIx-1)/2;

        if(!(dIy&1)) // even number
            dIy /= 2;
        else
            dIy = (dIy-1)/2;

        CVector<short,2> grad = { dIx, dIy };

        return grad;

}

CRGBImage::CRGBImage(size_t w, size_t h, unsigned char* data):
    CDenseArray<rgb>::CDenseArray(h,w) {

    unsigned char* pdata = reinterpret_cast<unsigned char*>(this->Data().get());

    memcpy(pdata,data,3*w*h);

}

#ifdef QT_GUI_LIB
CRGBImage::CRGBImage(const QImage& img):
    CDenseArray<rgb>::CDenseArray(img.height(),img.width()) {

    for(size_t i=0; i<NRows(); i++) {

        for(size_t j=0; j<NCols(); j++) {

            QRgb val = img.pixel(j,i);

            this->operator ()(i,j) = { (unsigned char)qRed(val),
                                       (unsigned char)qGreen(val),
                                       (unsigned char)qBlue(val) };

        }

    }


}


CRGBImage::operator QImage() const {

    QImage qimg(NCols(),NRows(),QImage::Format_RGB888);

    for(size_t i=0; i<NRows(); i++) {

        for(size_t j=0; j<NCols(); j++) {

            rgb val = CDenseArray<rgb>::Get(i,j);

            qimg.setPixel(j,i,qRgb(val.Get(0),val.Get(1),val.Get(2)));

        }

    }

    return qimg;

}
#endif

/*CRGBImage CRGBImage::Clone(CRectangle<int> roi) {

    CRGBImage result(roi.Width(),roi.Height());

    CVector<size_t,2> tl, br;
    tl = roi.TopLeft();
    br = roi.BottomRight();

    for(int i=tl.Get(1); i<=br.Get(1); i++) {

        for(int j=tl.Get(0); j<=br.Get(0); j++) {

            if(i>=0 && i<Height() && j>=0 && j<Width())
                result(i-tl.Get(1),j-tl.Get(0)) = CDenseArray<rgb>::Get(i,j);

        }

    }

    return result;

}


template<typename T>
CRGBImage CRGBImage::Clone(CRectangle<T> roi) {


    // round width and height
    size_t w, h;
    w = (size_t)(fabs(roi.Width())+0.5);
    h = (size_t)(fabs(roi.Height())+0.5);

    CRGBImage result(width,height);


    CVector<T,2> tl, br;
    tl = roi.TopLeft();
    br = roi.BottomRight();



    // PROBLEM here is the signed-ness of the width/height of the rectangle!
    // TODO: fix this!

}*/

} // end of namespace

