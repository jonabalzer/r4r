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

#ifndef R4RFEATURE_H_
#define R4RFEATURE_H_

#include <stdio.h>
#include <list>
#include <iostream>
#include <fstream>
#include <memory>

#include "rect.h"
#include "darray.h"
#include "vecn.h"


using namespace std;

namespace R4R {

// some forward declarations
class CAbstractDescriptor;
template<class Array> class CDescriptorAggregator;
template<class Array> class CSubsampleAggregator;

/*! \brief interest point in \f$\mathbb{R}^n\f$
 *
 *
 */
template<typename T,u_int n>
class CInterestPoint {

    friend class CAbstractDescriptor;
    template<class Array> friend class CDescriptorAggregator;
    template<class Array> friend class CSubsampleAggregator;

public:

    //! Standard constructor.
    CInterestPoint();

    //! Constructor.
    CInterestPoint(CVector<T,n>& location, float scale, T quality);

    //! Destructor.
    ~CInterestPoint();

    //! Returns the scale.
    float GetScale() const { return m_scale; }

    //! Returns the strength/quality of the feature.
    T GetQuality() const { return m_quality; }

    //! Sets the strength/quality of the feature.
    void SetQuality(T quality) { m_quality = quality; }

    //! Checks two features for equality.
    bool operator==(const CInterestPoint<T,n>& x);

    //! Checks if two features are not the same.
    bool operator!=(const CInterestPoint<T,n>& x) { return !operator==(x); }

    //! Adds a descriptor to the feature.
    void AttachDescriptor(const char* id, shared_ptr<CAbstractDescriptor> descriptor);

    //! Returns the number of descriptors attached to the feature.
    size_t NoDescriptors() const { return m_descriptors.size(); }

    //! Checks whether a descriptor of a given name exists.
    bool HasDescriptor(const char* name) const;

    //! Returns feature location w.r.t. to inherent scale.
    const CVector<T,n>& GetLocation() const { return m_location; }

    //! Returns feature location w.r.t. to native scale.
    CVector<T,n> GetLocationAtNativeScale() const;

    //! Access to the descriptor container.
    std::map<string,shared_ptr<CAbstractDescriptor> >& GetDescriptors() { return m_descriptors; }

    //! Looks for descriptor with a specified name. \TODO: Better return iterator.
    shared_ptr<CAbstractDescriptor> GetDescriptor(const char* name);

    //! Looks for a descriptor at specified position.
    shared_ptr<CAbstractDescriptor> GetDescriptor(u_int no);

    //! Looks for name of descriptor at specified position.
    string GetDescriptorName(u_int no);

    //! Writes a feature to an output stream.
    template<typename U,u_int m> friend std::ostream& operator <<(std::ostream& os, const CInterestPoint<U,m>& x);

    //! Writes a feature to a file stream.
    template<typename U,u_int m> friend std::ofstream& operator <<(std::ofstream& os, const CInterestPoint<U,m>& x);

    //! Reads a from a file stream.
    template<typename U,u_int m> friend std::ifstream& operator >>(std::ifstream& is, CInterestPoint<U,m>& x);

    /*! \brief Writes a set of features from file.
     * \param[in] filename file name
     * \param[in] features list of features
     * \param[in] comment comment
     * \returns true if the operation was successful
     */
    template<template<class U, class Allocator = std::allocator<U> > class Container>
    static bool SaveToFile(const char* filename, const Container<CInterestPoint<T,n> >& features, const char* comment = nullptr);

    /*! \brief Reads a set of features from file.
     * \param[in] filename file name
     * \param[out] features vector of features
     * \param[in] type floating point precision of the container
     * \returns number of features read, \f$-1\f$ on error
     */
    static int LoadFromFile(const char* filename, std::vector<CInterestPoint<T,n> >& features, string& comment);

private:

    CVector<T,n> m_location;                                                //! feature location
    float m_scale;                                                          //!< scale
    T m_quality;                                                            //!< quality measure, e.g., for storing motion estimation residual
    std::map<string,shared_ptr<CAbstractDescriptor> > m_descriptors;        //!< attached descriptors

};

} // end of namespace

#endif /* FEATURE_H_ */
