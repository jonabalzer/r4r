/*
 * tparams.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: jbalzer
 */

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

string CParameters::GetStringParameter(const char* name) {

	string param;

	try {

		if(m_string_params.find(name)==m_string_params.end())
			throw name;
		else
			param = m_string_params[name];

	}
	catch(const char* str) {

		cout << "ERROR: Parameter " << str << " was not loaded. Enter it now: ";

		cin >> param;

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

		cout << "ERROR: Parameter " << str << " was not loaded. Enter it now: ";

		cin >> param;

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

		cout << "ERROR: Parameter " << str << " was not loaded. Enter it now: ";



		cin >> param;

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

			}

		 	 break;

		 case 2:

		 	 {

		 	 int val;
		 	 in >> val;

		 	 m_int_params.insert(pair<string,int>(name,val));

		 	 }

		 	 break;

		 case 3:

		 	 {

			 double val;
		 	 in >> val;

		 	 m_double_params.insert(pair<string,double>(name,val));

		 	 }

		 	 break;

		 default:

			 break;

		 }

	 }

	in.close();

	//cout << (*this);

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

	return os;

}

}
