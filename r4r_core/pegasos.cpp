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

#include "pegasos.h"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>

using namespace std;

namespace R4R {

template <class Vector>
CPegasos<Vector>::CPegasos(vector<Vector>* features, vector<int>* labels, Vector w, double lambda, size_t k, size_t I):
	m_features(features),
	m_labels(labels),
	m_lambda(lambda),
	m_k(k),
	m_I(I),
	m_i(0),
	m_stepsize(1.0/m_lambda),
	m_w(w),
	m_b(0)
{

	m_res.push_back(CalcResidual(m_w,m_b));

}

template <class Vector>
double CPegasos<Vector>::CalcResidual(Vector w, double b) {

	double ws = 0.5*m_lambda*Kernel(w,w,b,b);

	double loss = 0;

	for(size_t i=0; i<m_features->size();i++)
		loss += max<double>(0,1.0-m_labels->at(i)*Kernel(w,m_features->at(i),b,1));

	loss /= (double)m_features->size();

	return ws + loss;

}

template <class Vector2>
ostream& operator << (ostream& os, CPegasos<Vector2>& x) {

	os << "-------------------------" << endl;
	os << "i = " << (int)x.m_i << endl;
	os << "Step size = " << x.m_stepsize << endl;
	os << "Residual = " << x.m_res[x.m_i] << endl;
	os << "w = " << endl << x.m_w << endl;
	os << "b = " << x.m_b << endl;
	os << "------------------------- " << endl;

	return os;
}

template <class Vector>
int CPegasos<Vector>::SubsampleFeatures(vector<Vector>* features, vector<int>* labels) {

	assert(m_k<m_features->size());

	for(size_t k=0; k<m_k; k++) {

        int rk = (int)(rand())%m_features->size();

		double margin = m_labels->at(rk)*Kernel(m_w,m_features->at(rk),m_b,1);

		if(margin<1) {

			features->push_back(m_features->at(rk));
			labels->push_back(m_labels->at(rk));

		}

	}

	return features->size();

}

template <class Vector>
bool CPegasos<Vector>::PerformDescentStep() {

	// ensure that the current solution is feasible
	assert(sqrt(Kernel(m_w,m_w,m_b,m_b))<1/sqrt(m_lambda));

	// check whether the maximum number of iterations has been reached
	if(m_i==m_I)
		return 1;

	// increment step variable
	m_i++;

	vector<Vector> features;
	vector<int> labels;

	// subsample features
	SubsampleFeatures(&features,&labels);

	// calculate gradient
	Vector grad(m_features->at(0));
	grad.Scale(0);
	double gradb = 0;

	for(size_t i=0;i<features.size();i++) {

		grad = grad + KernelGradW(m_w,features.at(i))*labels.at(i);
		gradb += labels.at(i)*KernelGradB(m_b,1);

	}

    grad.Scale(-1/(double)m_k);
	gradb /= -(double)m_k;

	grad = grad + m_w*m_lambda;
	gradb += m_lambda*m_b;

	// set step size
	SetArmijoStepSize(grad,gradb);

	// descent step
	grad.Scale(m_stepsize);
	m_w = m_w - grad;
	m_b += -m_stepsize*gradb;

	// projection step
	double s = min<double>(1.0,1/(sqrt(Kernel(m_w,m_w,m_b,m_b))*sqrt(m_lambda)));

	m_w.Scale(s);
	m_b *= s;

	// calc residual
    m_res.push_back(CalcResidual(m_w,m_b));

	return 0;

}

template <class Vector>
void CPegasos<Vector>::SetStepSize() {

	m_stepsize = 1/(m_lambda*((double)(m_i)));

}

template <class Vector>
void CPegasos<Vector>::Train(size_t n, bool silent) {

	// print initial state
	if(!silent)
		cout << *this << endl;

	for(size_t i=0;i<n;i++) {

		if(PerformDescentStep())
			break;

		if(!silent)
			cout << *this << endl;

	}

}

template <class Vector>
void CPegasos<Vector>::SetArmijoStepSize(Vector grad, double gradb) {

	SetStepSize();		// try learning rate first

	while(m_stepsize>0) {

		double flin, f;

		flin = CalcResidual(m_w,m_b) - (ARMIJO_C*m_stepsize)*Kernel(grad,grad,gradb,gradb);
		f = CalcResidual(m_w - grad*m_stepsize,m_b - m_stepsize*gradb);

		if(f<=flin)
			break;
		else
			m_stepsize *= ARMIJO_RHO;

	}

}

template <class Vector>
bool CPegasos<Vector>::SaveResidual(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file.\n";
		return 1;

	 }

	for(size_t k=0;k<m_res.size();k++)
		out << (int)k << " " << m_res[k] << endl;

	out.close();

	return 0;

}

template class CPegasos<vec>;
template class CPegasos<mat>;

template ostream& operator<< (ostream& os, CPegasos<vec>& x);
template ostream& operator<< (ostream& os, CPegasos<mat>& x);
template ostream& operator<< (ostream& os, CPegasos<smat>& x);

template <class Vector>
CPegasosMI<Vector>::CPegasosMI(vector<Vector>* features, vector<int>* labels, vector<int>* bags, Vector w, double lambda, size_t k, size_t I)
	:CPegasos<Vector>(features,labels,w,lambda,k,I),
	 m_bags(bags)
	 {


}

template <class Vector>
bool CPegasosMI<Vector>::OptimizePositiveBagHeuristic() {

	for(size_t i=0; i<m_labels->size();i++) {

		// check whether we are in a positive bag
		if(m_bags->at(i)>0) {

			// switch a label
			m_labels->at(i) *= -1;

			// calc the residual
			double res = CPegasos<Vector>::CalcResidual(m_w,m_b);

			// if bigger than actual value, switch back
			if(res>=m_res[m_i])
				m_labels->at(i) *= -1;

		}

	}

	return 0;

}

template <class Vector>
void CPegasosMI<Vector>::Train(size_t n, bool silent) {

	// print initial state
	if(!silent)
		cout << *this << endl;

	for(size_t i=0;i<n;i++) {

		if(CPegasos<Vector>::PerformDescentStep())
			break;

		OptimizePositiveBagHeuristic();

		if(!silent)
			cout << *this << endl;

	}

}

template class CPegasosMI<vec>;
template class CPegasosMI<mat>;


}
