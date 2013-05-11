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
