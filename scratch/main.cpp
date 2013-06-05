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


using namespace std;
using namespace R4R;
using namespace cv;

#include <opennurbs/opennurbs.h>

int main()
{




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



    CSplineCurve<double> curve = CSplineCurve<double>(2,3,5);
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
    ob4.close();
//    for(size_t i=0; i<cv.NCols(); i++) {

//        vec col = cv.GetColumn(i);
//        col.Scale(i);
//        cv.SetColumn(i,col);
//        //col(1)=0;

//    }

//    vec x0 = curve.Evaluate(1);

//    x0.Scale(3);
//   // x0(0)=3.4;
//   // x0(1)=1.2;

//    cout << "x0 " << x0 << endl;


//    CMercerKernel<double> kernel(2);
//    double t0, t1;
//    t0 = omp_get_wtime();
//    vec xc = curve.FindLocallyClosestPoint(x0,kernel,0.9,1e-4);

//    cout << xc << endl;

//    t1 = omp_get_wtime();
//    cout << t1-t0 << " s" << endl;
//    cout << xc << endl;


//   std::vector<CDescriptorFileHeader> headers;

//   ETYPE type;
//   cout << type << endl;

//   void* data;

//   data = CFeature::LoadDescriptors("/home/jbalzer/edkfsk;f.txt",headers,type);
//  // size_t stride = headers[0].NElems();


//   for(size_t k=0; k<headers.size(); k++) {

//        Mat img(headers[k].NCols(),headers[k].NRows(),CV_8UC1);

//        float* pdata = data+k*stride;

//        for(size_t i=0;i<headers[k].NRows();i++) {

//            for(size_t j=0;j<headers[k].NCols();j++) {

//                img.at<unsigned char>(j,i) = pdata[i*headers[k].NCols()+j];

//            }

//        }

//        imshow("Test",img);
//        waitKey();

//   }

   //ofstream mout("img.txt");
   //mout << img << endl;
   //mout.close();


//    CRectangle<double> roi = CRectangle<double>(0,0,1,1,0);
//    mat B(3,3);
//    B.Rand();
//    CDescriptor<mat>* desc = new CDescriptor<mat>(B);

//    cout << B << endl;

//    CFeature x(0,0,0,0);


//    shared_ptr<CDescriptor<mat> > ptr = shared_ptr<CDescriptor<mat> >(desc);
//    x.AttachDescriptor("ID",ptr);
//    x.AttachDescriptor("ID2",ptr);


//    list<CFeature> out, in;
//    out.push_back(x);
//    //out.push_back(x);

//    CFeature::SaveToFile("newtest.txt",out);

//    CFeature::OpenFromFile("newtest.txt",in);

//    list<CFeature>::iterator it;

//    for(it=in.begin(); it!=in.end(); it++) {

//        cout << *it << endl;

//        if(it->HasDescriptor("ID")) {

//            shared_ptr<CAbstractDescriptor> pdesc = it->GetDescriptor("ID");

//            cout << pdesc->GetType() << endl;
//            CDescriptor<mat>* desc = (CDescriptor<mat>*)pdesc.get();

//            mat& container = desc->Get();
//            cout << container << endl;


//        }

//        if(it->HasDescriptor("ID2")) {

//            shared_ptr<CAbstractDescriptor> pdesc = it->GetDescriptor("ID2");

//            cout << pdesc->GetType() << endl;
//            CDescriptor<mat>* desc = (CDescriptor<mat>*)pdesc.get();

//            mat& container = desc->Get();
//            cout << container << endl;


//        }


//    }




}

