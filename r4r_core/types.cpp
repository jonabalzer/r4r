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

#include "types.h"
#include "vecn.h"

namespace R4R {

template<typename T> ETYPE GetEType() { return ETYPE::NA; }
template<> ETYPE GetEType<bool>() { return ETYPE::B1U; }
template<> ETYPE GetEType<unsigned char>() { return ETYPE::C1U; }
template<> ETYPE GetEType<char>() { return ETYPE::C1S; }
template<> ETYPE GetEType<rgb>() { return ETYPE::C1U3; }
template<> ETYPE GetEType<unsigned short>() { return ETYPE::S2U; }
template<> ETYPE GetEType<short>() { return ETYPE::S2S; }
template<> ETYPE GetEType<int>() { return ETYPE::I4S; }
template<> ETYPE GetEType<unsigned int>() { return ETYPE::I4U; }
template<> ETYPE GetEType<float>() { return ETYPE::F4S; }
template<> ETYPE GetEType<long>() { return ETYPE::L8S; }
template<> ETYPE GetEType<unsigned long>() { return ETYPE::L8U; }
template<> ETYPE GetEType<double>() { return ETYPE::D8S; }
template<> ETYPE GetEType<vec3>() { return ETYPE::D8S3; }
template<> ETYPE GetEType<vec3f>() { return ETYPE::F4S3; }
template<> ETYPE GetEType<std::string>() { return ETYPE::STRING; }

// instantiations
template ETYPE GetEType<bool>();
template ETYPE GetEType<unsigned char>();
template ETYPE GetEType<char>();
template ETYPE GetEType<rgb>();
template ETYPE GetEType<unsigned short>();
template ETYPE GetEType<short>();
template ETYPE GetEType<int>();
template ETYPE GetEType<unsigned int>();
template ETYPE GetEType<float>();
template ETYPE GetEType<long>();
template ETYPE GetEType<unsigned long>();
template ETYPE GetEType<double>();
template ETYPE GetEType<vec3>();
template ETYPE GetEType<vec3f>();
template ETYPE GetEType<std::string>();



}
