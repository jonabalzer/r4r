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
#include <iostream>
#include <fstream>
#include <string>

#include "cam.h"
#include "rutils.h"
#include "factor.h"

using namespace std;

namespace R4R {

vector<vec2> CAbstractCam::Project(const vector<vec3>& x) const {

    vector<vec2> result(x.size());

#pragma omp parallel for
    for(size_t i=0; i<x.size(); i++)
        result[i] = Project(x[i]);

    return result;

}

vector<vec2f> CAbstractCam::Project(const vector<vec3f>& x) const {

    vector<vec2f> result(x.size());

#pragma omp parallel for
    for(size_t i=0; i<x.size(); i++)
        result[i] = Project(x[i]);

    return result;

}

vector<vec3> CAbstractCam::Normalize(const vector<vec2>& u) const {

    vector<vec3> result(u.size());

#pragma omp parallel for
    for(size_t i=0; i<u.size(); i++)
        result[i] = Normalize(u[i]);

    return result;

}

vector<vec3f> CAbstractCam::Normalize(const std::vector<vec2f> &u) const {

    vector<vec3f> result(u.size());

#pragma omp parallel for
    for(size_t i=0; i<u.size(); i++)
        result[i] = Normalize(u[i]);

    return result;

}

vector<vec2> CAbstractCam::Flow(const std::vector<vec3>& x, const std::vector<vec3>& dx) const {

    vector<vec2> result(x.size());

#pragma omp parallel for
    for(size_t i=0; i<x.size(); i++)
        result[i] = Flow(x[i],dx[i]);

    return result;

}

vector<vec2f> CAbstractCam::Flow(const std::vector<vec3f>& x, const std::vector<vec3f>& dx) const {

    vector<vec2f> result(x.size());

#pragma omp parallel for
    for(size_t i=0; i<x.size(); i++)
        result[i] = Flow(x[i],dx[i]);

    return result;

}

CPinholeCam::CPinholeCam() {

    fill_n(m_size,2,0);
    fill_n(m_f,2,0);
    fill_n(m_c,2,0);
    fill_n(m_k,5,0);
    m_alpha = 0;

}

CPinholeCam::CPinholeCam(double fu, double fv, double cu, double cv) {

    fill_n(m_size,2,0);
    m_f[0] = fu;
    m_f[1] = fv;
    m_c[0] = cu;
    m_c[1] = cv;
    fill_n(m_k,5,0);
    m_alpha = 0;

}

CPinholeCam::CPinholeCam(size_t w, size_t h, double fu, double fv, double cu, double cv) {

    m_size[0] = w;
    m_size[1] = h;
    m_f[0] = fu;
    m_f[1] = fv;
    m_c[0] = cu;
    m_c[1] = cv;
    fill_n(m_k,5,0);
    m_alpha = 0;

}


CPinholeCam::CPinholeCam(size_t w, size_t h) {

    m_size[0] = w;
    m_size[1] = h;
    m_f[0] = min(m_size[0],m_size[1]);
    m_f[1] = m_f[0];
    m_c[0] = 0.5*w;
    m_c[1] = 0.5*h;

    for(size_t i=0;i<5;i++)
        m_k[i]=0;

    m_alpha = 0;

}

vec2 CPinholeCam::Project(const vec3& x) const {

    vec2 xn = { x.Get(0)/x.Get(2), x.Get(1)/x.Get(2) };

    // radial distortion
    double r = xn.Norm2();

    vec2 dx;
    dx(0) = 2*m_k[2]*xn.Get(0)*xn.Get(1) + m_k[3]*(r*r + 2*xn.Get(0)*xn.Get(0));
    dx(1) = m_k[2]*(r*r + 2*xn.Get(1)*xn.Get(1)) + 2*m_k[3]*xn.Get(0)*xn.Get(1);

    double fac = 1 + m_k[0]*r*r + m_k[1]*r*r*r*r + m_k[4]*r*r*r*r*r*r;
    vec2 xd = xn*fac + dx;

    // transform to pixel coordinates
    vec2 xp;
    xp(0) = m_f[0]*(xd.Get(0) + m_alpha*xd.Get(1)) + m_c[0];
    xp(1) = m_f[1]*xd.Get(1) + m_c[1];

    return xp;

}

vec2f CPinholeCam::Project(const vec3f& x) const {

    vec2f xn = { x.Get(0)/x.Get(2), x.Get(1)/x.Get(2) };

    // radial distortion
    double r = xn.Norm2();

    vec2f dx;
    dx(0) = 2*m_k[2]*xn.Get(0)*xn.Get(1) + m_k[3]*(r*r + 2*xn.Get(0)*xn.Get(0));
    dx(1) = m_k[2]*(r*r + 2*xn.Get(1)*xn.Get(1)) + 2*m_k[3]*xn.Get(0)*xn.Get(1);

    double fac = 1 + m_k[0]*r*r + m_k[1]*r*r*r*r + m_k[4]*r*r*r*r*r*r;
    vec2f xd = xn*fac + dx;

    // transform to pixel coordinates
    vec2f xp;
    xp(0) = m_f[0]*(xd.Get(0) + m_alpha*xd.Get(1)) + m_c[0];
    xp(1) = m_f[1]*xd.Get(1) + m_c[1];

    return xp;

}


void CPinholeCam::Project(const vec3& x, vec2& u, mat& J) const {

    u = Project(x);

    J(0,0) = m_f[0]/x.Get(2);
    J(0,2) = -(u.Get(0)-m_c[0])/x.Get(2);
    J(1,1) = m_f[1]/x.Get(2);
    J(1,2) = -(u.Get(1)-m_c[1])/x.Get(2);

}

void CPinholeCam::Project(const vec3f& x, vec2f& u, matf& J) const {

    u = Project(x);

    J(0,0) = m_f[0]/x.Get(2);
    J(0,2) = -(u.Get(0)-m_c[0])/x.Get(2);
    J(1,1) = m_f[1]/x.Get(2);
    J(1,2) = -(u.Get(1)-m_c[1])/x.Get(2);

}

vec2 CPinholeCam::Flow(const vec3& x, const vec3& dx) const {

    vec2 u = Project(x);

    double J00 = m_f[0]/x.Get(2);
    double J02 = -(u.Get(0)-m_c[0])/x.Get(2);
    double J11 = m_f[1]/x.Get(2);
    double J12 = -(u(1)-m_c[1])/x.Get(2);

    vec2 du = { J00*dx.Get(0) + J02*dx.Get(2), J11*dx.Get(1) + J12*dx.Get(2) };

    return du;

}

vec2f CPinholeCam::Flow(const vec3f& x, const vec3f& dx) const {

    vec2f u = Project(x);

    float J00 = m_f[0]/x.Get(2);
    float J02 = -(u.Get(0)-m_c[0])/x.Get(2);
    float J11 = m_f[1]/x.Get(2);
    float J12 = -(u(1)-m_c[1])/x.Get(2);

    vec2f du = { J00*dx.Get(0) + J02*dx.Get(2), J11*dx.Get(1) + J12*dx.Get(2) };

    return du;

}


vec3 CPinholeCam::Normalize(const vec2& u) const {

    // FIXME: add undistortion
    vec3 x;
    x(0) = (1.0/m_f[0])*(u.Get(0)-m_c[0]);
    x(1) = (1.0/m_f[1])*(u.Get(1)-m_c[1]);
    x(2) = 1.0;

    return x;

}

vec3f CPinholeCam::Normalize(const vec2f& u) const {

    // TODO: add undistortion
    vec3f x;
    x(0) = (1.0/m_f[0])*(u.Get(0)-m_c[0]);
    x(1) = (1.0/m_f[1])*(u.Get(1)-m_c[1]);
    x(2) = 1.0;

    return x;

}

ostream& operator<< (ostream& os, const CPinholeCam& x) {

    os << "# dims" << endl;
    os << x.m_size[0] << " " << x.m_size[1] << endl;
    os << "# focal lengths" << endl;
    os << x.m_f[0] << " " << x.m_f[1] << endl;
    os << "# principle point" << endl;
    os << x.m_c[0] << " " << x.m_c[1] << endl;
    os << "# radial distortion coefficients" << endl;
    os << x.m_k[0] << " " << x.m_k[1] << " " << x.m_k[2] << " " << x.m_k[3] << " " << x.m_k[4] << endl;
    os << "# skew coefficient" << endl;
    os << x.m_alpha;

    return os;

}

istream& operator >> (istream& is, CPinholeCam& x) {

    string linebuffer;

    getline(is,linebuffer);

    is >> x.m_size[0];
    is >> x.m_size[1];
    is.get();

    getline(is,linebuffer);

    is >> x.m_f[0];
    is >> x.m_f[1];
    is.get();

    getline(is,linebuffer);

    is >> x.m_c[0];
    is >> x.m_c[1];
    is.get();

    getline(is,linebuffer);

    is >> x.m_k[0];
    is >> x.m_k[1];
    is >> x.m_k[2];
    is >> x.m_k[3];
    is >> x.m_k[4];
    is.get();

    getline(is,linebuffer);

    is >> x.m_alpha;

    return is;

}

bool CPinholeCam::operator==(CAbstractCam& cam) {

    CPinholeCam* pcam = dynamic_cast<CPinholeCam*>(&cam);

    bool result = (m_size[0]==pcam->m_size[0] &&
                   m_size[1]==pcam->m_size[1] &&
                   m_f[0]==pcam->m_f[0] &&
                   m_f[1]==pcam->m_f[1] &&
                   m_c[0]==pcam->m_c[0] &&
                   m_c[1]==pcam->m_c[1] &&
                   m_k[0]==pcam->m_k[0] &&
                   m_k[1]==pcam->m_k[1] &&
                   m_k[2]==pcam->m_k[2] &&
                   m_k[3]==pcam->m_k[3] &&
                   m_k[4]==pcam->m_k[4]);

    return result;

}

mat CPinholeCam::GetProjectionMatrix() const {

    mat P(3,3);

    P(0,0) = m_f[0];
    P(1,1) = m_f[1];
    P(0,2) = m_c[0];
    P(1,2) = m_c[1];
    P(2,2) = 1;

    return P;

}


template<typename T>
CView<T>::CView(CAbstractCam& cam, u_int index):
    m_cam(cam),
    m_F(),
    m_Finv(),
    m_index(index) {}

template<typename T>
CView<T>::CView(CAbstractCam& cam, const CRigidMotion<T,3> &F, u_int index):
    m_cam(cam),
    m_F(F),
    m_Finv(F),
    m_index(index) {

    m_Finv.Invert();

}

template<typename T>
CView<T>::CView(const CView<T>& view):
    m_cam(view.m_cam),
    m_F(view.m_F),
    m_Finv(view.m_Finv),
    m_index(view.m_index) {

}

template<typename T>
CView<T> CView<T>::operator=(const CView<T>& view) {

    if(!(m_cam==view.m_cam && m_F==view.m_F && m_Finv==view.m_Finv)) {

        m_cam = view.m_cam;
        m_F = view.m_F;
        m_Finv = view.m_Finv;
        m_index = view.m_index;

    }

    return *this;

}

template<typename T>
CVector<T,2> CView<T>::Project(const CVector<T,3>& x) const {

    CVector<T,3> xc = m_F.Transform(x);

    return m_cam.Project(xc);

}

template<typename T>
void CView<T>::Project(const CVector<T,3>& x, CVector<T,2>& u, CDenseArray<T>& J) const {

    CVector<T,3> xc = m_F.Transform(x);

    m_cam.Project(xc,u,J);

}

template<typename T>
CVector<T,2> CView<T>::Flow(const CVector<T,3>& x, const CVector<T,3>& dx) const {

    CVector<T,3> dxc = m_F.DifferentialTransform(dx);
    CVector<T,3> xc = m_F.Transform(x);

    return m_cam.Flow(xc,dxc);

}

template<typename T>
vector<CVector<T,2> > CView<T>::Project(const std::vector<CVector<T,3> >& x) const {

    vector<CVector<T,3> > xc = m_F.Transform(x);

    return m_cam.Project(xc);

}

template<typename T>
CVector<T,3> CView<T>::Normalize(const CVector<T,2>& u) const {

    CVector<T,3> xc = m_cam.Normalize(u);

    return m_Finv.Transform(xc);

}

template<typename T>
vector<CVector<T,3> > CView<T>::Normalize(const std::vector<CVector<T,2> >& u) const {

    vector<CVector<T,3> > xc = m_cam.Normalize(u);

    return m_Finv.Transform(xc);

}

template<typename T>
vector<CVector<T,2> > CView<T>::Flow(const vector<CVector<T,3> >& x, const vector<CVector<T,3> >& dx) const {

     vector<CVector<T,3> > dxc = m_F.DifferentialTransform(dx);
     vector<CVector<T,3> > xc = m_F.Transform(x);

     return m_cam.Flow(xc,dxc);

}

template<typename T>
bool CView<T>::SaveToFile(const char* filename) {

    ofstream out(filename);

    if(!out) {

        cout << "ERROR: Could not open file " << filename << "." << endl;
        return 1;

    }

    m_cam.Write(out);
    out << endl;
    out << m_F;

    out.close();

    return 0;

}

template<typename T>
bool CView<T>::OpenFromFile(const char* filename) {

    ifstream in(filename);

    if(!in) {

        cout << "ERROR: Could not open " << filename << "." << endl;
        return 1;

     }

    m_cam.Read(in);

    in.get();

    in >> m_F;

    in.close();

    m_Finv = m_F;
    m_Finv.Invert();

    return 0;

}

template<typename U>
ostream& operator << (ostream& os, const CView<U>& x) {

    // guarantees polymorphic behavior
    x.m_cam.Write(os);
    os << endl;

    os << x.m_F << endl;
    os << x.m_Finv;

    return os;

}

template <typename T>
CVector<T,3> CView<T>::GetLocation() const {

    // F is always world->cam, but the origin is translation of cam->world Finv
    return m_Finv.GetTranslation();

}

template <typename T>
CVector<T,3> CView<T>::GetPrincipalAxis() const {

    // F is always world->cam, but the origin is translation of cam->world Finv
    return m_Finv.GetPrincipalAxis();

}

template <typename T>
void CView<T>::Translate(const CVector<T,3>& t) {

    for(u_int i=0; i<3; i++)
        m_Finv(i,3) += t.Get(i);

    m_F = m_Finv;
    m_F.Invert();

}

template <typename T>
void CView<T>::DifferentialTranslate(const CVector<T,3>& t) {

    // first convert into world coordinates
    CVector<T,3> tw = m_Finv.DifferentialTransform(t);

    // now modify frames
    for(u_int i=0; i<3; i++)
        m_Finv(i,3) += tw.Get(i);

    m_F = m_Finv;
    m_F.Invert();

}

template <typename T>
void CView<T>::Orbit(const CVector<T,3>& center, const CVector<T,3>& axis) {

    // first get location in world
    CVector<T,3> origin = this->GetLocation();

    // lever
    CVector<T,3> d = origin - center;

    // create rotation
    CRigidMotion<T,3> R(0,0,0,axis.Get(0),axis.Get(1),axis.Get(2));

    // rotate lever and convert it back into point
    d = R.Transform(d);
    origin = center + d;

    // rotate axes
    CTransformation<T,3> Finv = R*m_Finv;
    m_Finv = reinterpret_cast<CRigidMotion<T,3>& >(Finv);

    // set new origin
    for(u_int i=0; i<3; i++)
        m_Finv(i,3) = origin.Get(i);

    // also set inverse
    m_F = m_Finv;
    m_F.Invert();

}



template class CView<float>;
template class CView<double>;
template ostream& operator << (ostream& os, const CView<float>& x);
template ostream& operator << (ostream& os, const CView<double>& x);


template <typename T>
bool CDistanceViewComparator<T>::operator ()(const CView<T>& x, const CView<T>& y) {

    CVector<T,3> locx = x.GetLocation();
    CVector<T,3> locy = y.GetLocation();

    return (locx - m_x).Norm2() < (locy - m_x).Norm2();

}

template class CDistanceViewComparator<double>;

template <typename T>
bool CAngleViewComparator<T>::operator ()(const CView<T>& x, const CView<T>& y) {

    // flip relation, cost is 1 - key
    return this->GetKey(x) > this->GetKey(y);

}

template <typename T>
T CAngleViewComparator<T>::GetKey(const CView<T>& x) {

    CVector<T,3> dx = x.GetLocation() - m_x;

    dx.Normalize();

    return InnerProduct(dx,m_n);

}

template class CAngleViewComparator<double>;

template <typename T>
bool CDistanceViewPairComparator<T>::operator ()(const std::pair<CView<T>,CView<T> >& x, const std::pair<CView<T>,CView<T> >& y) {

    CVector<T,3> dx = x.first.GetLocation() - x.second.GetLocation();
    CVector<T,3> dy = y.first.GetLocation() - y.second.GetLocation();

    return dx.Norm2() < dy.Norm2();

}

template class CDistanceViewPairComparator<double>;

template <typename T>
bool CAngleViewPairComparator<T>::operator ()(const std::pair<CView<T>,CView<T> >& x, const std::pair<CView<T>,CView<T> >& y) {

    // the principal axes are normalized!
    T ax = InnerProduct(x.first.GetPrincipalAxis(),x.second.GetPrincipalAxis());
    T ay = InnerProduct(y.first.GetPrincipalAxis(),y.second.GetPrincipalAxis());

    return ay > ax;  // same as 1 - ax < 1 - ay

}

template class CAngleViewPairComparator<double>;

}
