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

#include "params.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

using namespace std;

namespace R4R {

void CParameters::Set(const char* name, const char* val) {

	map<string,string>::iterator it = m_string_params.find(name);

	if(it!=m_string_params.end())
		m_string_params.erase(it);

	m_string_params.insert(pair<string,string>(name,val));

}

void CParameters::Set(const char* name, int val) {

	map<string,int>::iterator it = m_int_params.find(name);

	if(it!=m_int_params.end())
		m_int_params.erase(it);

	m_int_params.insert(pair<string,int>(name,val));

}

void CParameters::Set(const char* name, double val) {

	map<string,double>::iterator it = m_double_params.find(name);

	if(it!=m_double_params.end())
		m_double_params.erase(it);

	m_double_params.insert(pair<string,double>(name,val));

}

void CParameters::Set(const char* name, const vector<int>& val) {

    map<string,vector<int> >::iterator it = m_ints_params.find(name);

    if(it!=m_ints_params.end())
        m_ints_params.erase(it);

    m_ints_params.insert(pair<string,vector<int> >(name,val));

}

string CParameters::GetStringParameter(const char* name) {

	string param;

	try {

		if(m_string_params.find(name)==m_string_params.end())
			throw name;
		else
			param = m_string_params[name];

	}
	catch(const char* str) {

        cout << "ERROR: Parameter " << str << " was not loaded." << endl; // Enter it now: ";

        //cin >> param;
        param = string("-1");

        Set(name,param.c_str());

	}

	return param;

}

int CParameters::GetIntParameter(const char* name) {

	int param;

	try {

		if(m_int_params.find(name)==m_int_params.end())
			throw name;
		else
			param = m_int_params[name];

	}
	catch(const char* str) {

        cout << "ERROR: Parameter " << str << " was not loaded." << endl; // Enter it now: ";

        //cin >> param;
        param = -1;

		Set(name,param);

	}

	return param;

}

double CParameters::GetDoubleParameter(const char* name) {

	double param;

	try {

		if(m_double_params.find(name)==m_double_params.end())
			throw name;
		else
			param = m_double_params[name];

	}
	catch(const char* str) {

        cout << "ERROR: Parameter " << str << " was not loaded." << endl; // Enter it now: ";

        //cin >> param;
        param = -1;

		Set(name,param);


	}

	return param;

}

vector<int> CParameters::GetIntsParameter(const char* name) {

    vector<int> param;

    try {

        if(m_ints_params.find(name)==m_ints_params.end())
            throw name;
        else
            param = m_ints_params[name];

    }
    catch(const char* str) {

        cout << "ERROR: Parameter " << str << " was not loaded." << endl;
        //cout << "Enter list of integers. Stop input by entering -1:" << endl;

        /*int input;

        while(true) {

            cin >> input;

            if(input==-1)
              break;
            else
              param.push_back(input);

        }*/

        param.push_back(-1);

        Set(name,param);

    }

    return param;

}

bool CParameters::SaveToFile(const char* filename) {

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

        case 4:

        {

            // read size
            size_t n;
            in >> n;

            vector<int> ints;

            for(size_t i=0; i<n; i++) {

                int val;
                in >> val;
                ints.push_back(val);

            }

            m_ints_params.insert(pair<string,vector<int> >(name,ints));

        }

        default:

            break;

        }

    }

    in.close();

    return 0;

}

ostream& operator<< (ostream& os, CParameters& x) {

	map<string,string>::iterator it;

	for(it=x.m_string_params.begin(); it!=x.m_string_params.end(); it++)
		os << it->first << " " << 1 << " " << it->second << endl;

	map<string,int>::iterator it2;

	for(it2=x.m_int_params.begin(); it2!=x.m_int_params.end(); it2++)
		os << it2->first << " " << 2 << " " << it2->second << endl;

	map<string,double>::iterator it3;

	for(it3=x.m_double_params.begin(); it3!=x.m_double_params.end(); it3++)
		os << it3->first << " " << 3 << " " << it3->second << endl;

    map<string,vector<int> >::iterator it4;

    for(it4=x.m_ints_params.begin(); it4!=x.m_ints_params.end(); it4++) {

        os << it4->first << " " << 4 << " " << it4->second.size();

        for(size_t i=0; i<it4->second.size(); i++)
            os << " " << it4->second.at(i);

        os << endl;

    }

	return os;

}

}
