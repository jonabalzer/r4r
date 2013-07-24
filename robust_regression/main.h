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

#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <math.h>


#include "lm.h"

using namespace std;
using namespace R4R;

class CRosenbrockFunction:public CLeastSquaresProblem<mat,double> {

public:

    //! Constructor.
    CRosenbrockFunction();

    //! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&)
    void ComputeResidualAndJacobian(vec& r, mat& J);

    //! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
    void ComputeResidual(vec& r);

    //! Prints out the ground truth parameters.
    vec GetGroundTruth() { vec x(2); x(0) = 1; x(1) = 1; return x; }

private:

};

class COsbourneFunction:public CLeastSquaresProblem<mat,double> {

public:

    //! Constructor.
    COsbourneFunction();

    //! \copydoc CLeastSquaresProblem::ComputeResidualAndJacobian(vec&,Matrix&)
    void ComputeResidualAndJacobian(vec& r, mat& J);

    //! \copydoc CLeastSquaresProblem::ComputeResidual(vec&)
    void ComputeResidual(vec& r);

    //! Generates abscissae and ordinates pairs.
    void DisturbSamplePoints(size_t noutlier, double strength);

    //! Prints out the data points.
    void PrintData();

    //! Prints out the ground truth parameters.
    vec GetGroundTruth();

    //! Evaluates the function itself.
    static double EvalFunction(double t);

    //! Save data to file.
    bool SaveData(const char* filename);

private:

    static double m_y[33];
    static double m_t[33];

};




#endif // MAIN_H
