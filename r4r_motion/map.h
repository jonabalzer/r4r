/*
 * map.h
 *
 *  Created on: Jan 30, 2013
 *      Author: jbalzer
 */

#ifndef MAP_H_
#define MAP_H_

#include "feature.h"

namespace R4R {


class CMap: public std::list<CFeature> {

public:

	//! Standard constructor.
	CMap();


	//! Saves map to a PLY file.
	bool SaveToFile(const char* filename);

private:


};


}

#endif /* MAP_H_ */
