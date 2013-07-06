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

#ifndef R4RPARAMS_H_
#define R4RPARAMS_H_

#include <stdlib.h>
#include <string>
#include <map>
#include <vector>

namespace R4R {


/*! \brief container for user-definable parameters
 *
 *
 *
 */
class CParameters {

public:

    //! Read parameters from file.
	bool OpenFromFile(const char* filename);

	//! Write parameters to file.
	bool SaveToFile(const char* filename);

    //! Writes all parameters to an output stream.
	friend std::ostream& operator << (std::ostream& os, CParameters& x);

	//! Sets a string parameter.
	void Set(const char* name, const char* val);

	//! Sets an integer parameter.
	void Set(const char* name, int val);

	//! Sets a float parameter.
    void Set(const char* name, double val);

    //! Sets an integer vector parameter.
    void Set(const char* name, const std::vector<int>& val);

	//! Gets a string parameter.
	std::string GetStringParameter(const char* name);

	//! Gets an integer parameter.
	int GetIntParameter(const char* name);

	//! Gets a float parameter.
	double GetDoubleParameter(const char* name);

    //! Gets an integer vector parameter.
    std::vector<int> GetIntsParameter(const char* name);

private:

    std::map<std::string,int> m_int_params;					//!< container for integer parameters
    std::map<std::string,double> m_double_params;			//!< container for real parameters
    std::map<std::string,std::string> m_string_params;		//!< container for string parameters
    std::map<std::string,std::vector<int> > m_ints_params;  //!< container for series of integers

};

}

#endif /* PARAMS_H_ */
