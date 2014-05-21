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

#ifndef R4RTYPES_H_
#define R4RTYPES_H_

namespace R4R {

/*! \brief
 * a list of data types commonly used in R4R
 *
 */
enum class ETYPE {  NA = 0,
                    B1U = 1,
                    C1U = 2, C1S = 3, C1U3 = 102,
                    S2U = 4, S2S = 5,
                    I4S = 6, I4U = 7,
                    F4S = 8, F4S3 = 108,
                    L8S = 9, L8U = 10,
                    D8S = 11, D8S3 = 111,
                    STRING = 19 };

template<typename T> ETYPE GetEType();

}

#endif /* TYPES_H_ */
