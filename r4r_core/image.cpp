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


template<typename T>
CIntImage<T>::CIntImage(size_t width, size_t height):
    CDenseArray<T>::CDenseArray(height,width) {}

template<typename T>
void CIntImage<T>::Compute() {

    if(this->NCols()>=this->NRows()) {

        for(size_t d=1; d<this->NRows(); d++) {

            // process the quadratic part of the matrix
            for(size_t i=0; i<d; i++) {

                if(i==0) {

                    // column above (d+d)
                    this->operator ()(d,i) += this->Get(d-1,i);

                    // row left of (d+d)
                    this->operator ()(i,d) += this->Get(i,d-1);

                }
                else {

                    // column above (d+d)
                    this->operator ()(d,i) += this->Get(d,i-1) + this->Get(d-1,i) - this->Get(d-1,i-1);

                    // row left of (d+d)
                    this->operator ()(i,d) += this->Get(i,d-1) + this->Get(i-1,d) - this->Get(i-1,d-1);

                }

            }

            // value at (d+d)
            this->operator ()(d,d) += this->Get(d,d-1) + this->Get(d-1,d) - this->Get(d-1,d-1);

        }

        // process the remaining columns
        for(size_t j=this->NRows(); j<this->NCols(); j++) {

            for(size_t i=0; i<this->NRows(); i++) {

                if(i==0)
                    this->operator ()(i,j) += this->Get(i,j-1);
                else
                    this->operator ()(i,j) += this->Get(i,j-1) + this->Get(i-1,j) - this->Get(i-1,j-1);

            }

        }

    }
    else {

        for(size_t d=1; d<this->NCols(); d++) {

            for(size_t i=0; i<d; i++) {

                if(i==0) {

                    // column above (d+d)
                    this->operator ()(d,i) += this->Get(d-1,i);

                    // row left of (d+d)
                    this->operator ()(i,d) += this->Get(i,d-1);

                }
                else {

                    // column above (d+d)
                    this->operator ()(d,i) += this->Get(d,i-1) + this->Get(d-1,i) - this->Get(d-1,i-1);

                    // row left of (d+d)
                    this->operator ()(i,d) += this->Get(i,d-1) + this->Get(i-1,d) - this->Get(i-1,d-1);

                }


            }

            // value at (d+d)
            this->operator ()(d,d) += this->Get(d,d-1) + this->Get(d-1,d) - this->Get(d-1,d-1);

        }

        // process the remaining row
        for(size_t i=this->NCols(); i<this->NRows(); i++) {

            for(size_t j=0; j<this->NCols(); j++) {

                if(j==0)
                    this->operator ()(i,j) += this->Get(i-1,j);
                else
                    this->operator ()(i,j) += this->Get(i,j-1) + this->Get(i-1,j) - this->Get(i-1,j-1);


            }

        }

    }

}

template<typename T>
void CIntImage<T>::AddMass(const CVector<double,2>& x, T val) {

    double dx, dy;
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

template<typename T>
void CIntImage<T>::AddMass(const CVector<size_t,2>& x, T val) {

    this->operator ()(x.Get(1),x.Get(0)) += val;

}

template<typename T>
template<typename U>
U CIntImage<T>::Evaluate(const CVector<double,2>& x, const CVector<double,2>& hsize) const {

    CVector<double,2> tl = this->ProjectToBoundary(x - hsize);
    CVector<double,2> br = this->ProjectToBoundary(x + hsize);
    CVector<double,2> bl = { tl.Get(0), br.Get(1) };
    CVector<double,2> tr = { br.Get(0), tl.Get(1) };

    return this->template Get<U>(tl) + this->template Get<U>(br) - this->template Get<U>(bl) - this->template Get<U>(tr);

}

template double CIntImage<double>::Evaluate<double>(const CVector<double,2>& x, const CVector<double,2>& hsize) const;

template<typename T>
T CIntImage<T>::EvaluateApproximately(const CVector<double,2>& x, const CVector<double,2>& hsize) const {

    CVector<double,2> tl = this->ProjectToBoundary(x - hsize);
    CVector<double,2> br = this->ProjectToBoundary(x + hsize);
    CVector<double,2> bl = { tl.Get(0), br.Get(1) };
    CVector<double,2> tr = { br.Get(0), tl.Get(1) };

    return this->Get(tl) + this->Get(br) - this->Get(bl) - this->Get(tr);

}

template class CIntImage<size_t>;
template class CIntImage<double>;
template class CIntImage<float>;


} // end of namespace

