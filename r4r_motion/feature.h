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

#ifndef R4RFEATURE_H_
#define R4RFEATURE_H_

#include <stdio.h>
#include <list>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>

#include "rect.h"
#include "darray.h"
#include "vecn.h"

namespace R4R {

// some forward declarations
class CAbstractDescriptor;

/*! \brief interest point in \f$\mathbb{R}^n\f$
 *
 *
 */
template<typename T,u_int n>
class CInterestPoint {

    friend class CAbstractDescriptor;

public:

    //! Standard constructor.
    CInterestPoint();

    //! Constructor.
    CInterestPoint(std::initializer_list<T> ilist);

    //! Constructor.
    CInterestPoint(const CVector<T,n>& location);

    //! Constructor.
    CInterestPoint(const CVector<T,n>& location, float scale, T quality);

    /*! \brief Copy constructor.
     *
     * \param[in] name name of descriptor to include in the copy
     *
     * All but one descriptor are ignored during construction.
     *
     */
    CInterestPoint(const CInterestPoint<T,n>& x, std::string name);

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
    void AttachDescriptor(const char* id, std::shared_ptr<CAbstractDescriptor> descriptor);

    //! Returns the number of descriptors attached to the feature.
    size_t NoDescriptors() const { return m_descriptors.size(); }

    //! Checks whether a descriptor of a given name exists.
    bool HasDescriptor(const char* name) const;

    //! Returns feature location w.r.t. to inherent scale.
    const CVector<T,n>& GetLocation() const { return m_location; }

    //! Sets location of the feature..
    void SetLocation(const CVector<T,n>& location) { m_location = location; }

    //! Returns feature location w.r.t. to native scale.
    CVector<T,n> GetLocationAtNativeScale() const;

    //! Access to the descriptor container.
    std::unordered_map<std::string,std::shared_ptr<CAbstractDescriptor> >& GetDescriptors() { return m_descriptors; }

    //! Looks for descriptor with a specified name. \TODO: Better return iterator.
    const std::shared_ptr<CAbstractDescriptor>& GetDescriptor(const char* name) const;

    //! Looks for a descriptor at specified position.
    const std::shared_ptr<CAbstractDescriptor>& GetDescriptor(u_int no) const;

    //! Looks for name of descriptor at specified position.
    std::string GetDescriptorName(u_int no);

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
    static int LoadFromFile(const char* filename, std::vector<CInterestPoint<T,n> >& features, std::string& comment);

private:

    CVector<T,n> m_location;                                                          //! feature location
    float m_scale;                                                                    //!< scale
    T m_quality;                                                                      //!< quality measure, e.g., for storing motion estimation residual
    std::unordered_map<std::string,std::shared_ptr<CAbstractDescriptor> > m_descriptors;        //!< attached descriptors

};

typedef CInterestPoint<float,2> imfeature;

//template<typename T>
//using C2dFeature = CInterestPoint<T,2>;

//template<typename T>
//using C3dFeature = CInterestPoint<T,3>;


} // end of namespace

#endif /* FEATURE_H_ */
