/*
 * map.cpp
 *
 *  Created on: Jan 30, 2013
 *      Author: jbalzer
 */



#include "map.h"
#include <fstream>

using namespace std;

namespace R4R {


CMap::CMap():list<CFeature>()  {}


bool CMap::SaveToFile(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file " << filename << "." << endl;
		return 1;

	 }

	out << "ply" << endl;
	out <<  "format ascii 1.0" << endl;
	out <<  "comment" << endl;
	out << "element vertex " << size() << endl;
	out << "property float32 x" << endl;
	out << "property float32 y" << endl;
	out << "property float32 z" << endl;
	out << "end_header" << endl;

	list<CFeature>::iterator it;

	for(it=begin(); it!=end(); it++) {

		vec pt = it->GetLocation();
		out << pt.Get(0) << " " << pt.Get(1) << " " << pt.Get(2) << endl;

	}

	out.close();

	return 0;

}

}
