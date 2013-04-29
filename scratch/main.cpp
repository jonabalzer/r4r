#include <iostream>
#include <fstream>
#include "darray.h"
#include "descriptor.h"
#include "feature.h"
#include <omp.h>
#include <stdio.h>
#include "splinecurve.h"

using namespace std;
using namespace R4R;
using namespace cv;
#include <opennurbs/opennurbs.h>

int main()
{



    ON::Begin();

    ON_TextLog dump;

    int order =3 ;
    ON_NurbsCurve curve2 = ON_NurbsCurve(3,0,order,13);

    curve2.MakeClampedUniformKnotVector(0.1);



    double N[order*order];

    for(size_t i=0; i<curve2.m_cv_count; i++) {

        ON_3dPoint pt;
        pt.x = i;
        pt.y = i;
        pt.z = 0;
        curve2.SetCV(i,pt);

    }
    curve2.Dump(dump);
ON_3dPoint pt = curve2.PointAt(0.05);

ON_3dVector tangent = curve2.TangentAt(0.05);

cout << pt.x << " " << pt.y << " " << pt.z << " " << endl;

cout << tangent.x << " " << tangent.y << " " << tangent.z << " " << endl;



   ON::End();


    CSplineCurve<double> curve = CSplineCurve<double>(2,2,13);

    curve.MakeClampedUniformKnotVector(0,1);
    mat& cv = curve.GetCVData();

    cv.Ones();


    curve.Print();

    for(size_t i=0; i<cv.NCols(); i++) {

        vec col = cv.GetColumn(i);
        col.Scale(i);
        //col(1)=0;

    }
        cout << cv << endl;
//    cout << cv << endl;

   mat x = curve.Tangent(0.95);


    //vec x = curve.Evaluate(0.05);
    cout << x << endl;







    //}
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

