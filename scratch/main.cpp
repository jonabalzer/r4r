#include <iostream>
#include <fstream>
#include "darray.h"
#include "descriptor.h"
#include "feature.h"
#include <omp.h>
#include <stdio.h>
#include "splinecurve.h"
#include "factor.h"
#include "types.h"
#include "trafo.h"
#include "cam.h"

using namespace std;
using namespace R4R;
using namespace cv;

bool rw_descriptors();

#include <opennurbs/opennurbs.h>

int main()
{

    mat M(3,4);
    M.Rand(0,1);

    cout << M << endl;
    CRigidMotion<double,3> F(M);

    mat R = F.GetJacobian();
    cout << R << endl;

    CVector<double,3> t=F.GetTranslation();
    cout << t << endl;

    vec t2(t);
    cout << "t2" << t2 << endl;

   /* CPinholeCam cam;
    CView myview(cam,F);
    myview.OpenFromFile("camtest.txt");


    cout << myview << endl;

    vec3 x={1,4,2};

    vec2 y = myview.Project(x);
    cout << y << endl;

    vector<vec3> in;
    for(size_t i=0; i<1000000; i++)
        in.push_back(x);

    double t0 = omp_get_wtime();
    vector<vec2> out = myview.Project(in);
    double t1 = omp_get_wtime();
    cout << t1-t0 << endl;

    cout << out[0] << endl;*/




    //rw_descriptors();
    //cout << "Returned" << endl;

//    ON::Begin();

//    ON_TextLog dump;

//    int ncp = 13;
//    int order = 4;
//    int p = order - 1;
//    int nk = ncp + p - 1;

//    cout << "Nk: " << nk << endl;
//    ON_NurbsCurve curve2 = ON_NurbsCurve(3,0,order,ncp);

//    double delta = 1.0/(nk - 2*(order-2) - 1);

//    curve2.MakeClampedUniformKnotVector(delta);

//   // double N[order*order];

//    for(size_t i=0; i<curve2.m_cv_count; i++) {

//        ON_3dPoint pt;
//        pt.x = i;
//        pt.y = i;
//        pt.z = 0;
//        curve2.SetCV(i,pt);

//    }

//    curve2.Dump(dump);

//    ON_3dPoint pt;// = curve2.PointAt(0.05);

//    ON_3dVector tangent;// = curve2.TangentAt(0.05);

//    //curve2.Ev1Der(0.05,pt,tangent);

//    double v[6];
//    curve2.Evaluate(0.05,1,3,&v[0]);
//    cout << v[0] << " " << v[1] << " " << v[2] << endl;
//    cout << v[3] << " " << v[4] << " " << v[5] << endl;

//    double N[order*order];
//    ON_EvaluateNurbsBasis(order,curve2.m_knot,0.05,N);
//    ON_EvaluateNurbsBasisDerivatives(order,curve2.m_knot,1, N);
//    cout << N[order+0] << " " << N[order+1] << " " << N[order+2] << " " << N[order+3] << endl;

//    int span = ON_NurbsSpanIndex(order,curve2.m_cv_count,curve2.m_knot,0.95,-1,0);

//    cout << "Span: " << span << endl;

//    //cout << pt.x << " " << pt.y << " " << pt.z << " " << endl;
//    //cout << tangent.x << " " << tangent.y << " " << tangent.z << " " << endl;



//    ON::End();



   /* CSplineCurve<double> curve = CSplineCurve<double>(2,3,5);
    curve.MakeClampedUniformKnotVector(0,1);
    mat& cv = curve.GetCVData();

    cv.Ones();

    double N[4*4];

    size_t n = 100;
    double dt = 1.0/(double(n-1));

    double b0[n];
    double b1[n];
    double b2[n];
    double b3[n];
    double b4[n];

    vec& knot = curve.GetKnotVector();
    curve.Print();
    for(size_t i=0; i<n; i++) {

        int span = curve.GetSpan(i*dt);



        CSplineCurve<double>::EvaluateNurbsBasis(4,knot.Data().get()+span,i*dt,N);

        cout << N[0] << " " << N[1] << " " << N[2] << " " << N[3] << endl;
        if(span==0) {

            b0[i]=N[0];
            b1[i]=N[1];
            b2[i]=N[2];
            b3[i]=N[3];
            b4[i] = 0;
        }
        else
        {



            b0[i]=0;
            b1[i]=N[0];
            b2[i]=N[1];
            b3[i]=N[2];
            b4[i]=N[3];



        }



    }


    ofstream ob0("b0.txt");
    ofstream ob1("b1.txt");
    ofstream ob2("b2.txt");
    ofstream ob3("b3.txt");
    ofstream ob4("b4.txt");

    for(size_t i=0; i<n; i++) {
        ob0 << b0[i] << " ";
        ob1 << b1[i] << " ";
        ob2 << b2[i] << " ";
        ob3 << b3[i] << " ";
        ob4 << b4[i] << " ";
    }

    ob0.close();
    ob1.close();
    ob2.close();
    ob3.close();
    ob4.close();*/





}


bool rw_descriptors() {

    // create interest point
    CVector<float,2> x = { 2.4, 2.0 };
    CInterestPoint<float,2> feature(x,3,2.0);
    cout << feature << endl;

    // create descriptor data, memory is allocated here
    mat B(3,3);
    B.Rand(0,1);
    cout << B << endl;

    // the container member of descriptor gets another reference to the data.
    CDescriptor<mat>* desc = new CDescriptor<mat>(B);

    // a shared pointer to the descriptor object is attached to feature
    shared_ptr<CDescriptor<mat> > ptr = shared_ptr<CDescriptor<mat> >(desc);
    feature.AttachDescriptor("ID",ptr);

    list<CInterestPoint<float,2> > out;
    out.push_back(feature);
    //out.push_back(feature);

    CInterestPoint<float,2>::SaveToFile("iotest.txt",out);

    vector<CInterestPoint<float,2> > in;

    string com;
    CInterestPoint<float,2>::LoadFromFile("iotest.txt",in,com);

    cout << com << endl;
    cout << in.size() << endl;

    for(size_t i=0; i<in.size(); i++) {

        cout << in[i] << endl;

        if(in[i].HasDescriptor("ID")) {

            shared_ptr<CAbstractDescriptor> pdesc = in[i].GetDescriptor("ID");

            cout << (int)pdesc->GetType() << endl;
            CDescriptor<mat>* desc = (CDescriptor<mat>*)pdesc.get();

            mat& container = desc->Get();
            cout << container << endl;

        }

    }

    return 0;

}


