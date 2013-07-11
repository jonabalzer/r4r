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

CRGBImage::CRGBImage(size_t w, size_t h, unsigned char* data):
    CDenseArray<rgb>::CDenseArray(h,w) {

    unsigned char* pdata = reinterpret_cast<unsigned char*>(this->Data().get());

    memcpy(pdata,data,3*w*h);

}

template<typename T>
CVector<T,3> CRGBImage::Get(const CVector<T,2>& p) {

    double u, v;
    u = p.Get(0);
    v = p.Get(1);

    int i = (int)floor(v);
    int j = (int)floor(u);

    if(i<0 || i>=NRows()-1 || j<0 || j>=NCols()-1)
        return vec3();

    double vd = v - i;
    double ud = u - j;

    vec3 A00, A01, A10, A11, A0, A1;
    A00 = vec3(CDenseArray<rgb>::Get(i,j));
    A01 = vec3(CDenseArray<rgb>::Get(i,j+1));
    A10 = vec3(CDenseArray<rgb>::Get(i+1,j));
    A11 = vec3(CDenseArray<rgb>::Get(i+1,j+1));

    A0 = A00*(1-ud) + A01*ud;
    A1 = A10*(1-ud) + A11*ud;

    return A0*(1-vd) + A1*vd;

}

template CVector<double,3> CRGBImage::Get<double>(const CVector<double,2>& p);
template CVector<float,3> CRGBImage::Get<float>(const CVector<float,2>& p);

template<typename T>
CVector<T,3> CRGBImage::Gradient(const CVector<T,2>& p, bool dir) {

    CVector<T,2> dp;

    if(dir)
        dp = { 0, 1 };
    else
        dp = { 1, 0 };

    CVector<T,3> I0, I1;
    I1 = Get(p+dp);
    I0 = Get(p-dp);

    return 0.5*(I1-I0);

}

template CVector<double,3> CRGBImage::Gradient<double>(const CVector<double,2>& p, bool dir);
template CVector<float,3> CRGBImage::Gradient<float>(const CVector<float,2>& p, bool dir);


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

}

