#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <math.h>


#include "lm.h"

using namespace std;
using namespace R4R;

class CRosenbrockFunction:public CLeastSquaresProblem<mat> {

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

class COsbourneFunction:public CLeastSquaresProblem<mat> {

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
