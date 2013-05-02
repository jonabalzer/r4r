//////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////////

#ifndef R4RIMAGE_H
#define R4RIMAGE_H

#include "types.h"

namespace R4R {

/*! \brief R4R's own gray value image class
 *
 *
 *
 */
class CImage: public CDenseArray<unsigned char> {

public:

    //! Constructor.
    CImage();

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

    //! Constructor.
    CRGBImage(size_t w, size_t h):CDenseArray<rgb>::CDenseArray(h,w){}

};

}

#endif // IMAGE_H
