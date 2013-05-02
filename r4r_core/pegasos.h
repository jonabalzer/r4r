/*
 * pegasos.h
 *
 *  Created on: Mar 14, 2012
 *      Author: jbalzer
 */

#ifndef R4RPEGASOS_H_
#define R4RPEGASOS_H_

#define ARMIJO_RHO 0.8
#define ARMIJO_C 0.25

#include <stdio.h>
#include "types.h"

namespace R4R {

/*! \brief PEGASOS stochastic subgradient SVM solver.
 *
 *
 *
 */
template<class Vector>
class CPegasos {

public:

	//! Constructor.
	CPegasos(std::vector<Vector>* features, std::vector<int>* labels, Vector w, double lambda = 1.0, size_t k = 1, size_t I = 100);

	//! Prints the state of all variables to an output stream.
	template<class Vector2> friend std::ostream& operator << (std::ostream& os, CPegasos<Vector2>& x);

	//! Triggers execution of a given number of descent steps.
	void Train(size_t n, bool silent = false);

	//! Calculates the objective function value.
	double CalcResidual(Vector w, double b);

	//! Saves the residual progression to a file.
	bool SaveResidual(const char* filename);

	//! Provides access to #m_w.
	Vector GetW() { return m_w; };

	//! Provides access to #m_b.
	double GetBias() { return m_b; }

	//! Classification with trained SVM.
	double Classify(Vector x) {	return Kernel(m_w,x,m_b,1); };

protected:

	// external variables
	std::vector<Vector>* m_features;			//!< pointer to a container of feature vectors
	std::vector<int>* m_labels;					//!< pointer to a container of labels
	double m_lambda;							//!< regularization parameter
	size_t m_k;									//!< size of mini patches
	size_t m_I;									//!< maximum number of iterations

	// internal states
	size_t m_i;									//!< iteration index
	double m_stepsize;							//!< step size
	Vector m_w;									//!< normal of separating hyperplane
	double m_b;									//!< bias term
	std::vector<double>	m_res;					//!< container for residual values over #m_i

	//! Randomly selects #m_k samples from all features.
	int SubsampleFeatures(std::vector<Vector>* features, std::vector<int>* labels);

	//! Sets step size according to learning rate.
	void SetStepSize();

	//! Armijo-type step size rule.
	void SetArmijoStepSize(Vector grad, double gradb);

	//! Executes a single descent step.
	bool PerformDescentStep();

	/*! \brief Evaluates kernel function.
	 *
	 * Overload this function to generalize the method to nonlinear kernels.
     *
     * \FIXME: No, it is more complicated than that.
	 *
	 */
    virtual double Kernel(Vector x, Vector y, double a, double b) { return Vector::InnerProduct(x,y) + a*b; };

	/*! Gradient of the kernel w.r.t. y in the subspace of #m_w.
	 *
	 * Overload this function to generalize the method to nonlinear kernels.
	 *
	 */
	virtual Vector KernelGradW(Vector x, Vector y) { return y; };

	/*! Partial derivative of kernel w.r.t. b in the subspace of #m_b.
	 *
	 * Overload this function to generalize the method to nonlinear kernels.
	 *
	 */
	virtual double KernelGradB(double a, double b) { return b; };

};

/*! \brief PEGASOS stochastic subgradient SVM-MIL solver
 *
 *
 *
 */
template<class Vector>
class CPegasosMI:public CPegasos<Vector> {

public:

	//! Constructor.
	CPegasosMI(std::vector<Vector>* features, std::vector<int>* labels, std::vector<int>* bags, Vector w, double lambda = 1.0, size_t k = 1, size_t I = 100);

	//! Invoke training.
	void Train(size_t n, bool silent = false);

private:

	std::vector<int>* m_bags;																//!< pointer to a container of bag indicators

	bool OptimizePositiveBagHeuristic();

	using CPegasos<Vector>::m_labels;
	using CPegasos<Vector>::m_w;
	using CPegasos<Vector>::m_b;
	using CPegasos<Vector>::m_res;
	using CPegasos<Vector>::m_i;


};


}

#endif /* PEGASOS_H_ */
