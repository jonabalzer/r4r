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

#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits>

#include "trafo.h"
#include "factor.h"

using namespace std;

namespace R4R {

template <typename T,u_int n>
CTransformation<T,n>::CTransformation() {

    fill_n(m_F,n*(n+1),0);

    for(u_int i=0; i<n; i++)
        m_F[n*i+i] = 1;

}

template <typename T,u_int n>
CTransformation<T,n>::CTransformation(const CDenseArray<T>& F) {

    assert(F.NRows()==n && F.NCols()==(n+1));

    memcpy(&m_F[0],F.Data().get(),n*(n+1)*sizeof(T));

}

template <typename T,u_int n>
CTransformation<T,n> CTransformation<T,n>::operator*(const CTransformation<T,n>& x) const {

    CTransformation<T,n> result;

    // iterate through cols of x
    for(u_int k=0; k<=n; k++) {

        for(u_int i=0; i<n; i++) {

            result.m_F[k*n+i] = 0;

            for(u_int j=0; j<n; j++)
                result.m_F[k*n+i] += m_F[n*j+i]*x.m_F[k*n+j];

            if(k==n)
                result.m_F[k*n+i] += m_F[n*n+i];

        }

    }

    return result;

}

template <typename T,u_int n>
T CTransformation<T,n>::Get(u_int i, u_int j) const {

    // col major
    return m_F[n*j + i];

}

template <typename T,u_int n>
T& CTransformation<T,n>::operator ()(u_int i, u_int j) {

    // col major
    return m_F[n*j + i];

}

template <typename T,u_int n>
CTransformation<T,n>::operator CDenseArray<T>() const {

    CDenseArray<T> result(n+1,n+1);
    result(n,n) = 1;

    for(u_int i=0; i<n; i++) {

        for(u_int j=0; j<=n; j++)
            result(i,j) = m_F[n*j + i];

    }

    return result;

}

template <typename T,u_int n>
CVector<T,n> CTransformation<T,n>::Transform(const CVector<T, n>& x) const {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++) {

        for(u_int j=0; j<n; j++)
            result(i) += m_F[n*j+i]*x.Get(j);

        result(i) += m_F[n*n+i];

    }

    return result;

}

template <typename T,u_int n>
CVector<T,n> CTransformation<T,n>::DifferentialTransform(const CVector<T,n>& x) const {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++) {

        for(u_int j=0; j<n; j++)
            result(i) += m_F[n*j+i]*x.Get(j);

    }

    return result;

}

template <typename T,u_int n>
CVector<T,n> CTransformation<T,n>::Transform(const T* x) {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++) {

        for(u_int j=0; j<n; j++)
            result(i) += m_F[n*j+i]*x[j];

    result(i) += m_F[n*n+i];

    }

    return result;

}

template <typename T,u_int n>
vector<CVector<T,n> > CTransformation<T,n>::Transform(const vector<CVector<T,n> >& x) const {

    vector<CVector<T,n> > result(x.size());

#pragma omp parallel for
    for(size_t i=0; i<x.size(); i++)
        result[i] = Transform(x[i]);

    return result;

}

template <typename T,u_int n>
vector<CVector<T,n> > CTransformation<T,n>::DifferentialTransform(const vector<CVector<T,n> >& x) {

    vector<CVector<T,n> > result(x.size());

#pragma omp parallel for
    for(size_t i=0; i<x.size(); i++)
        result[i] = DifferentialTransform(x[i]);

    return result;

}


template <class U,u_int m>
ostream& operator << (ostream& os, const CTransformation<U,m>& x) {

    os << "# coordinate transformation" << endl;

    for(u_int i=0; i<m; i++) {

        for(u_int j=0; j<m; j++)
            os << x.m_F[m*j+i] << " ";

        if(i<m-1)
            os << x.m_F[m*m+i] << endl;
        else
            os << x.m_F[m*m+i];

    }

    return os;

}


template <class U,u_int m>
istream& operator >> (istream& is, CTransformation<U,m>& x) {

    // get comment line
    string linebuffer;
    getline(is,linebuffer);

    // read transformation
    for(u_int i=0; i<m; i++) {

        for(u_int j=0; j<(m+1); j++)
            is >> x.m_F[m*j+i];

    }

    return is;

}

template <typename T,u_int n>
bool CTransformation<T,n>::Invert() {

    T* temp = new T[n*(n+1)];
    memcpy(temp,m_F,n*(n+1)*sizeof(T));

    if(n==3) {

        // inverse of det
        T det = temp[0]*(temp[4]*temp[8] - temp[5]*temp[7])
                - temp[3]*(temp[8]*temp[1] - temp[2]*temp[7])
                + temp[6]*(temp[1]*temp[5] - temp[2]*temp[4]);

        if(det==0)
            return 1;

        T deti = 1/det;

        // first linear part
        m_F[0] = deti*(temp[4]*temp[8] - temp[5]*temp[7]);
        m_F[1] = -deti*(temp[1]*temp[8] - temp[2]*temp[7]);
        m_F[2] = deti*(temp[1]*temp[5] - temp[2]*temp[4]);
        m_F[3] = -deti*(temp[3]*temp[8] - temp[5]*temp[6]);
        m_F[4] = deti*(temp[0]*temp[8] - temp[2]*temp[6]);
        m_F[5] = -deti*(temp[0]*temp[5] - temp[2]*temp[3]);
        m_F[6] = deti*(temp[3]*temp[7] - temp[4]*temp[6]);
        m_F[7] = -deti*(temp[0]*temp[7] - temp[1]*temp[6]);
        m_F[8] = deti*(temp[0]*temp[4] - temp[1]*temp[3]);

    }
    else {

        T det = temp[0]*temp[3]-temp[1]*temp[2];

        if(det==0)
            return 1;

        T deti = 1/det;

        // first linear part
        m_F[0] = deti*temp[3];
        m_F[1] = -deti*temp[1];
        m_F[2] = -deti*temp[2];
        m_F[3] = deti*temp[0];

    }

    // then affine part
    for(u_int i=0; i<n; i++) {

        m_F[n*n+i] = 0;

        for(u_int j=0; j<n; j++)
            m_F[n*n+i] -= m_F[n*j+i]*temp[n*n+j];

    }

    delete [] temp;

    return 0;

}

template <typename T,u_int n>
CDenseArray<T> CTransformation<T,n>::GetJacobian() {

    CDenseArray<T> result(n,n);

    memcpy(result.Data().get(),&m_F[0],n*n*sizeof(T));

    return result;

}

template <typename T,u_int n>
CVector<T,n> CTransformation<T,n>::GetTranslation() {

    CVector<T,n> result;

    for(u_int i=0; i<n; i++)
        result(i) = m_F[n*n+i];

    return result;

}

template <typename T,u_int n>
bool CTransformation<T,n>::operator ==(const CTransformation<T,n>& x) {

    bool result = true;

    for(size_t i=0; i<n*(n+1); i++)
        result *= (m_F[i]==x.m_F[i]);

    return result;

}


template class CTransformation<double,2>;
template class CTransformation<float,2>;
template class CTransformation<double,3>;
template class CTransformation<float,3>;
template ostream& operator << (ostream& os, const CTransformation<double,3>& x);
template ostream& operator << (ostream& os, const CTransformation<double,2>& x);
template ostream& operator << (ostream& os, const CTransformation<float,3>& x);
template ostream& operator << (ostream& os, const CTransformation<float,2>& x);
template istream& operator >> (istream& is, CTransformation<double,3>& x);
template istream& operator >> (istream& is, CTransformation<double,2>& x);
template istream& operator >> (istream& is, CTransformation<float,3>& x);
template istream& operator >> (istream& is, CTransformation<float,2>& x);

template<typename T>
CRotation<T,2>::CRotation(T o) {

    m_F[0] = cos(o);
    m_F[1] = sin(o);
    m_F[2] = -sin(o);
    m_F[3] = cos(o);

}

template<typename T>
bool CRotation<T,2>::Invert() {

    // transpose R
    swap(m_F[1],m_F[2]);

    return 0;

}

template class CRotation<double,2>;
template class CRotation<float,2>;

template<typename T>
CDifferentialRotation<T,2>::CDifferentialRotation(T o) {

    m_F[0] = -sin(o);
    m_F[1] = cos(o);
    m_F[2] = -cos(o);
    m_F[3] = -sin(o);

}

template class CDifferentialRotation<double,2>;
template class CDifferentialRotation<float,2>;

template<typename T>
CRotation<T,3>::CRotation(T o1, T o2, T o3) {

     Rodrigues(o1,o2,o3,m_F);

}

template<typename T>
void CRotation<T,3>::Rodrigues(const  T& o1, const T& o2, const T& o3, T* R) {

    T theta = sqrt(o1*o1+o2*o2+o3*o3);

    if(theta<numeric_limits<T>::epsilon()) {

        fill_n(R,9,0);

        R[0] = 1;
        R[4] = 1;
        R[8] = 1;

        return;

    }

    T ox, oy, oz;
    ox = o1/theta;
    oy = o2/theta;
    oz = o3/theta;

    T oxox, oxoy, oxoz, oyoy, oyoz, ozoz;
    oxox = ox*ox;
    oxoy = ox*oy;
    oxoz = ox*oz;
    oyoy = oy*oy;
    oyoz = oy*oz;
    ozoz = oz*oz;

    T sth, cth, mcth;
    sth = sin(theta);
    cth  = cos(theta);
    mcth = 1 - cth;

    R[0] = 1 - mcth*(oyoy+ozoz);    R[3] = -sth*oz + mcth*oxoy;		R[6] = sth*oy + mcth*oxoz;
    R[1] = sth*oz + mcth*oxoy;      R[4] = 1 - mcth*(ozoz + oxox);	R[7] = -sth*ox + mcth*oyoz;
    R[2] = - sth*oy + mcth*oxoz;	R[5] = sth*ox + mcth*oyoz;    	R[8] = 1 - mcth*(oxox+oyoy);

}

template<typename T>
void CRotation<T,3>::Log(const T* R, T& o1, T& o2, T& o3) {

    T arg = 0.5*((R[0]+R[4]+R[8])-1.0);

    if(fabs(arg)>1-numeric_limits<T>::epsilon())
        return;

    T theta = acos(arg);
    T s = theta/(2*sin(theta));

    o1 = -s*(R[7]-R[5]);
    o2 = s*(R[6]-R[2]);
    o3 = -s*(R[3]-R[1]);

}

template<typename T>
bool CRotation<T,3>::Invert() {

    // transpose R
    swap(m_F[1],m_F[3]);
    swap(m_F[2],m_F[6]);
    swap(m_F[5],m_F[7]);

    return 0;

}

template class CRotation<double,3>;
template class CRotation<float,3>;

template <typename T>
CDifferentialRotation<T,3>::CDifferentialRotation(T o1, T o2, T o3, u_int dim) {

    switch(dim) {

    case 0:
        Rodrigues(o1,o2,o3,nullptr,m_F,nullptr,nullptr);
        break;
    case 1:
        Rodrigues(o1,o2,o3,nullptr,nullptr,m_F,nullptr);
        break;
    case 2:
        Rodrigues(o1,o2,o3,nullptr,nullptr,nullptr,m_F);
        break;

    }

}

template <typename T>
void CDifferentialRotation<T,3>::Rodrigues(const T& o1, const T& o2, const T& o3, T* R, T* DRo1, T* DRo2, T* DRo3) {

    T theta = sqrt(o1*o1+o2*o2+o3*o3);

    if(theta<numeric_limits<T>::epsilon()) {

        if(R!=nullptr) {

            fill_n(R,9,0);
            R[0] = 1;
            R[4] = 1;
            R[8] = 1;

        }

        if(DRo1!=nullptr) {

            fill_n(DRo1,9,0);
            DRo1[5] = 1;
            DRo1[7] = -1;

        }

        if(DRo2!=nullptr) {

            fill_n(DRo2,9,0);
            DRo2[2] = -1;
            DRo2[6] = 1;

        }

        if(DRo3!=nullptr) {

            fill_n(DRo3,9,0);
            DRo3[1] = 1;
            DRo3[3] = -1;

        }

        return;

    }

    T ox, oy, oz;
    ox = o1/theta;
    oy = o2/theta;
    oz = o3/theta;

    T oxox, oxoy, oxoz, oyoy, oyoz, ozoz;
    oxox = ox*ox;
    oxoy = ox*oy;
    oxoz = ox*oz;
    oyoy = oy*oy;
    oyoz = oy*oz;
    ozoz = oz*oz;

    T sth, cth, mcth;
    sth = sin(theta);
    cth  = cos(theta);
    mcth = 1 - cth;

    if(R!=nullptr) {

        R[0] = 1 - mcth*(oyoy+ozoz);    R[3] = -sth*oz + mcth*oxoy;		R[6] = sth*oy + mcth*oxoz;
        R[1] = sth*oz + mcth*oxoy;      R[4] = 1 - mcth*(ozoz + oxox);	R[7] = -sth*ox + mcth*oyoz;
        R[2] = - sth*oy + mcth*oxoz;	R[5] = sth*ox + mcth*oyoz;    	R[8] = 1 - mcth*(oxox+oyoy);

    }

    T a, b, c, d;
    a =  sth/theta;
    b = mcth/theta;
    c = cth - a;
    d = sth - 2*b;

    if(DRo1!=nullptr) {

        DRo1[0] = -d*(oyoy + ozoz)*ox;			DRo1[3] = b*oy - c*oxoz + d*oxoy*ox;	DRo1[6] = b*oz + c*oxoy + d*oxoz*ox;
        DRo1[1] = b*oy + c*oxoz + d*oxoy*ox;	DRo1[4] = -2*b*ox - d*(ozoz + oxox)*ox; DRo1[7] = -a - c*oxox + d*oyoz*ox;
        DRo1[2] = b*oz - c*oxoy + d*oxoz*ox; 	DRo1[5] = a + c*oxox + d*oyoz*ox;		DRo1[8] = -2*b*ox - d*(oyoy + oxox)*ox;

    }

    if(DRo2!=nullptr) {

        DRo2[0] = -2*b*oy - d*(oyoy + ozoz)*oy;	DRo2[3] = b*ox - c*oyoz + d*oxoy*oy;	DRo2[6] = a + c*oyoy + d*oxoz*oy;
        DRo2[1] = b*ox + c*oyoz + d*oxoy*oy;	DRo2[4] = -d*(ozoz + oxox)*oy;			DRo2[7] =  b*oz - c*oxoy + d*oyoz*oy;
        DRo2[2] = -a - c*oyoy + d*oxoz*oy;		DRo2[5] =  b*oz + c*oxoy + d*oyoz*oy;	DRo2[8] = -2*b*oy - d*(oyoy + oxox)*oy;

    }

    if(DRo3!=nullptr) {

        DRo3[0] = -2*b*oz - d*(oyoy + ozoz)*oz;	DRo3[3] = -a - c*ozoz + d*oxoy*oz;		DRo3[6] = b*ox + c*oyoz + d*oxoz*oz;
        DRo3[1] = a + c*ozoz + d*oxoy*oz;		DRo3[4] = -2*b*oz - d*(ozoz + oxox)*oz; DRo3[7] = b*oy - c*oxoz + d*oyoz*oz;
        DRo3[2] =  b*ox - c*oyoz + d*oxoz*oz;	DRo3[5] = b*oy + c*oxoz + d*oyoz*oz;	DRo3[8] = - d*(oyoy + oxox)*oz;

    }

}

template class CDifferentialRotation<double,3>;
template class CDifferentialRotation<float,3>;

template<typename T>
CRigidMotion<T,2>::CRigidMotion(T t1, T t2, T o) {

    m_F[0] = cos(o);
    m_F[1] = sin(o);
    m_F[2] = -sin(o);
    m_F[3] = cos(o);
    m_F[4] = t1;
    m_F[5] = t2;

}


template<typename T>
bool CRigidMotion<T,2>::Invert() {

    // transpose R
    swap(m_F[1],m_F[2]);

    // -R'*t
    T temp[2];
    temp[0] = m_F[4];
    temp[1] = m_F[5];

    m_F[4] = -(m_F[0]*temp[0] + m_F[2]*temp[1]);
    m_F[5] = -(m_F[1]*temp[0] + m_F[3]*temp[1]);

    return 0;

}

template class CRigidMotion<double,2>;
template class CRigidMotion<float,2>;

template<typename T>
CRigidMotion<T,3>::CRigidMotion(T t1, T t2, T t3, T o1, T o2, T o3) {

    CRotation<T,3>::Rodrigues(o1,o2,o3,m_F);
    m_F[9] = t1;
    m_F[10] = t2;
    m_F[11] = t3;

}

template<typename T>
CRigidMotion<T,3>::CRigidMotion(const CDenseVector<T>& m) {

    // since constructor delegation is not yet supported by gcc 4.6
    m_F[9] = m.Get(0);
    m_F[10] = m.Get(1);
    m_F[11] = m.Get(2);
    CRotation<T,3>::Rodrigues(m.Get(3),m.Get(4),m.Get(5),m_F);


}

template<typename T>
bool CRigidMotion<T,3>::Invert() {

    // transpose R
    swap(m_F[1],m_F[3]);
    swap(m_F[2],m_F[6]);
    swap(m_F[5],m_F[7]);

    // -R'*t
    T temp[3];
    temp[0] = m_F[9];
    temp[1] = m_F[10];
    temp[2] = m_F[11];

    m_F[9] = -(m_F[0]*temp[0] + m_F[3]*temp[1] + m_F[6]*temp[2]);
    m_F[10] = -(m_F[1]*temp[0] + m_F[4]*temp[1] + m_F[7]*temp[2]);
    m_F[11] = -(m_F[2]*temp[0] + m_F[5]*temp[1] + m_F[8]*temp[2]);

    return 0;

}


template class CRigidMotion<double,3>;
template class CRigidMotion<float,3>;

}
