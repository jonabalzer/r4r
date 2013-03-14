#include <iostream>
#include <fstream>
#include "darray.h"
#include "feature.h"
#include "basic.h"
#include "rect.h"
#include "tracklet.h"
#include "params.h"
#include "descspecial.h"

using namespace std;
using namespace R4R;
using namespace cv;

int main()
{


    Mat img = imread("/home/jbalzer/Data/imgs/grador0.png");


    float di = 255.0/img.cols;

    for(size_t i=0; i<img.rows; i++ ) {

        for(size_t j=0; j<img.cols; j++) {


            img.at<unsigned char>(i,j) = (unsigned char)(j*di);

        }

    }

    vec u(2);
    u(0)=320;
    u(1)=240;

    namedWindow("Test",1);
    imshow("Test",img);
    waitKey(0);

    CRectangle<double> roi(u(0),u(1),8,8);
    //CFourierModulusDescriptor desc(roi,25);

    CHistogramOfGradients desc(roi);

    desc.Compute(img);

    vecf& test = desc.Get();


    cout << test << endl;
    //ofstream out("fft.txt");
    //out << test;
    //out.close();

//    CRectangle<double> roi = CRectangle<double>(0,0,1,1,0);

//    CIdentityDescriptor* desc = new CIdentityDescriptor(roi,0,7);

//    matf& B = desc->Get();

//    B(1,1)=20;





//    //cout << test << endl;

//    shared_ptr<CIdentityDescriptor> ptr = shared_ptr<CIdentityDescriptor>(desc);

//    x.AttachDescriptor("ID",ptr);

//    //x.AttachDescriptor("BRIEF",ptr);
//    //x.AttachDescriptor("ID",ptr);

//    list<CFeature> l;

//    CFeature y(2,3,0,0);
//    y.AttachDescriptor("ID",ptr);

//    l.push_back(x);
//    l.push_back(y);

//    CFeature::SaveDescriptors("bw.txt",l,"ID","Comment2");
//    vector<CDescriptorFileHeader> featin;

//    float* buffer = new float[10000000];

//    CFeature::LoadDescriptors("bw.txt",featin,buffer);

//    for(size_t i=0;i<featin.size(); i++)
//        cout << featin[i] << endl;


//    cout << buffer[0] << endl;

//    matf A(15,15,buffer);
//    cout << A << endl;


    //    CFeature::OpenFromFileBlockwise("bw.txt",featin,names,data,(15*15));

//    for(size_t i=0; i<names.size(); i++) {
//        cout << featin[i] << endl;
//        cout << names[i] << endl;

//    }

    //cout << data << endl;

    return 0;


}

