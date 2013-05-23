/*
 * trafo.h
 *
 *  Created on: Oct 16, 2012
 *      Author: jbalzer
 */

#ifndef TRAFO_H_
#define TRAFO_H_

#include "types.h"

namespace R4R {

/*! \brief transformation
 *
 * \todo Introduce template parameter for data type.
 *
 */
template<size_t size>
class CTransformation {

public:

	//! Constructor.
	CTransformation();

	//! Constructor.
	CTransformation(const mat& F);

	//! Invert transformation.
	virtual void Invert() {};

	//! Returns upper-left part of the homogeneous matrix.
	mat GetLinearPart() const;

	//! Returns translation vector.
	vec GetTranslation() const;

	//! Sets upper-left part of the homogeneous matrix.
	void SetLinearPart(const mat& A);

	//! Sets translation vector.
	void SetTranslation(const vec& t);

protected:

	mat m_F;

};



/*! \brief rotation
 *
 *
 *
 */
template<size_t size>
class CRotation:public CTransformation<size> {

public:

	//! Standard constructor.
	CRotation():CTransformation<size>() {};

	/*! \brief Constructor.
	 *
	 * \details Takes a matrix and projects it to \f$\mathrm{SO}(3)\f$ by Procrustes analysis.
	 *
	 */
	CRotation(const mat& A);

	//! Invert transformation.
	virtual void Invert() {};

protected:

	using CTransformation<size>::m_F;


};


template<> class CRotation<2>:public CTransformation<2> {

public:

	//! Constructor.
	CRotation(double o);

	//! Rodrigues formula.
	static mat Rodrigues(double o);

protected:

	using CTransformation<2>::m_F;

};

template<> class CRotation<3>:public CTransformation<3> {

public:

	//! Constructor.
	CRotation(double o1, double o2, double o3);

	//! Constructor.
	CRotation(vec omega);

	//! Rodrigues formula.
	static mat Rodrigues(double o1, double o2, double o3);

	//! Logarithm.
	static vec Log(const mat& R);

	//! Logarithm.
	vec Log();

protected:

	using CTransformation<3>::m_F;


};


/*! \brief rigid motion
 *
 *
 *
 */
template<size_t size>
class CRigidMotion:public CTransformation<size> {

public:

	//! Constructor.
	CRigidMotion():CTransformation<size>() {};

};





}



#endif /* TRAFO_H_ */
