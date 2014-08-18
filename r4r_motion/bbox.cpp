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

#include "bbox.h"

using namespace std;

namespace R4R {

template<typename T>
CBoundingBox<T>::CBoundingBox(const CVector<T,3>& lower, const CVector<T,3>& upper):
    m_lower(lower),
    m_upper(upper) {}


template<typename T>
vector<CVector<T,3> > CBoundingBox<T>::Corners() const {

    vector<CVector<T,3> > result(8);

    result[0] = m_lower;
    result[1] = { m_upper.Get(0), m_lower.Get(1), m_lower.Get(2) };
    result[2] = { m_lower.Get(0), m_upper.Get(1), m_lower.Get(2) };;
    result[3] = { m_upper.Get(0), m_upper.Get(1), m_lower.Get(2) };;
    result[4] = m_upper;
    result[5] = { m_upper.Get(0), m_lower.Get(1), m_upper.Get(2) };
    result[6] = { m_lower.Get(0), m_upper.Get(1), m_upper.Get(2) };;
    result[7] = { m_upper.Get(0), m_upper.Get(1), m_upper.Get(2) };;

    return result;

}

template class CBoundingBox<double>;
template class CBoundingBox<float>;

}
