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

#ifndef R4RBBOX_H
#define R4RBBOX_H

#include <vector>

#include "vecn.h"


namespace R4R {

template<typename T>
class CBoundingBox {

public:

    //! Constructor.
    CBoundingBox():m_lower(), m_upper() {}

    //! Constructor.
    CBoundingBox(const R4R::CVector<T,3>& lower, const R4R::CVector<T,3>& upper);

    //! Get corners.
    std::vector<R4R::CVector<T,3> > Corners() const;

    //! Access to the lower corner.
    const CVector<T,3>& Lower() { return m_lower; }

    //! Access to the upper corner.
    const CVector<T,3>& Upper() { return m_upper; }

    //! Barycenter.
    CVector<T,3>  Barycenter() { return (m_upper + m_lower)*0.5; }

private:

    R4R::CVector<T,3> m_lower;          //!< lower corner
    R4R::CVector<T,3> m_upper;          //!< upper corner

};

}

#endif // BBOX_H
