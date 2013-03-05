/*
 * tparams.h
 *
 *  Created on: Mar 26, 2012
 *      Author: jbalzer
 */

#ifndef PARAMS_H_
#define PARAMS_H_

#include <stdlib.h>
#include <string>
#include <map>

namespace R4R {


/*! \brief container for user-definable parameters
 *
 *
 *
 */
class CParameters {

public:

	std::map<std::string,int> m_int_params;					//!< container for integer parameters
	std::map<std::string,double> m_double_params;			//!< container for real parameters
	std::map<std::string,std::string> m_string_params;		//!< container for string parameters

	//! Read parameters from file.
	bool OpenFromFile(const char* filename);

	//! Write parameters to file.
	bool SaveToFile(const char* filename);

	// ! Writes all parameters to an output stream.
	friend std::ostream& operator << (std::ostream& os, CParameters& x);

	//! Sets a string parameter.
	void Set(const char* name, const char* val);

	//! Sets an integer parameter.
	void Set(const char* name, int val);

	//! Sets a float parameter.
	void Set(const char* name, double val);

	//! Gets a string parameter.
	std::string GetStringParameter(const char* name);

	//! Gets an integer parameter.
	int GetIntParameter(const char* name);

	//! Gets a float parameter.
	double GetDoubleParameter(const char* name);

};

}

#endif /* PARAMS_H_ */
