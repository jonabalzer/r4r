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

#include "params.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

using namespace std;

namespace R4R {

void CParameters::Set(const char* name, const char* val) {

    unordered_map<string,string>::iterator it = m_string_params.find(name);

	if(it!=m_string_params.end())
		m_string_params.erase(it);

	m_string_params.insert(pair<string,string>(name,val));

}

void CParameters::Set(const char* name, int val) {

    unordered_map<string,int>::iterator it = m_int_params.find(name);

	if(it!=m_int_params.end())
		m_int_params.erase(it);

	m_int_params.insert(pair<string,int>(name,val));

}

void CParameters::Set(const char* name, double val) {

    unordered_map<string,double>::iterator it = m_double_params.find(name);

	if(it!=m_double_params.end())
		m_double_params.erase(it);

	m_double_params.insert(pair<string,double>(name,val));

}

string CParameters::GetStringParameter(const char* name) const {

	string param;

	try {

        unordered_map<string,string>::const_iterator it = m_string_params.find(name);

        if(it==m_string_params.end())
			throw name;
		else
            param = it->second;

	}
	catch(const char* str) {

        cerr << "ERROR: Parameter " << str << " was not loaded." << endl;
        param = string("");

	}

	return param;

}

int CParameters::GetIntParameter(const char* name) const {

	int param;

	try {

        unordered_map<string,int>::const_iterator it = m_int_params.find(name);

        if(it==m_int_params.end())
			throw name;
		else
            param = it->second;

	}
	catch(const char* str) {

        cerr << "ERROR: Parameter " << str << " was not loaded." << endl;
        param = -1;

	}

	return param;

}

double CParameters::GetDoubleParameter(const char* name) const {

	double param;

	try {

        unordered_map<string,double>::const_iterator it = m_double_params.find(name);

        if(it==m_double_params.end())
			throw name;
		else
            param = it->second;

	}
	catch(const char* str) {

        cout << "ERROR: Parameter " << str << " was not loaded." << endl;
        param = -1.0;

	}

	return param;

}

bool CParameters::SaveToFile(const char* filename) const {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file.\n";
		return 1;

	 }

	out << *this;

	out.close();

	return 0;

}

bool  CParameters::OpenFromFile(const char* filename) {

    ifstream in(filename);

    if(!in) {

        cout << "ERROR: Could not open " << filename << "." << endl;
        return 1;

    }

    while (in.good()) {

        string name;
        int type;

        in >> name;
        in >> type;

        switch (type) {

        case 1:

        {

            string val;
            in >> val;

            m_string_params.insert(pair<string,string>(name,val));

            break;

        }

        case 2:

        {

            int val;
            in >> val;

            m_int_params.insert(pair<string,int>(name,val));

            break;

        }

        case 3:

        {

            double val;
            in >> val;

            m_double_params.insert(pair<string,double>(name,val));

            break;

        }

        default:

            break;

        }

    }

    in.close();

    return 0;

}

ostream& operator<< (ostream& os, const CParameters& x) {

    unordered_map<string,string>::const_iterator it;

	for(it=x.m_string_params.begin(); it!=x.m_string_params.end(); it++)
		os << it->first << " " << 1 << " " << it->second << endl;

    unordered_map<string,int>::const_iterator it2;

	for(it2=x.m_int_params.begin(); it2!=x.m_int_params.end(); it2++)
		os << it2->first << " " << 2 << " " << it2->second << endl;

    unordered_map<string,double>::const_iterator it3;

	for(it3=x.m_double_params.begin(); it3!=x.m_double_params.end(); it3++)
		os << it3->first << " " << 3 << " " << it3->second << endl;

	return os;

}

}
