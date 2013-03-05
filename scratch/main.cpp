#include <iostream>
#include <fstream>
#include "darray.h"
#include "feature.h"
#include "basic.h"
#include "rect.h"
#include "tracklet.h"
#include "params.h"

using namespace std;
using namespace R4R;


int main()
{


    vec u(2);
    u(0)=1;
    u(1)=3;
    CFeature x = CFeature(u,3,0);

    CRectangle<double> roi = CRectangle<double>(0,0,1,1,0);

    //CBRIEF* desc = new CBRIEF(roi);

    CIdentityDescriptor* desc = new CIdentityDescriptor(roi,0,7);


    //cout << test << endl;

    shared_ptr<CIdentityDescriptor> ptr = shared_ptr<CIdentityDescriptor>(desc);

    x.AttachDescriptor("ID",ptr);

    //x.AttachDescriptor("BRIEF",ptr);
    //x.AttachDescriptor("ID",ptr);

    list<CFeature> l;

    CFeature y(2,3,0,0);
    y.AttachDescriptor("ID",ptr);

    l.push_back(x);
    l.push_back(y);

    CFeature::SaveToFileBlockwise("bw.txt",l,(15*15),0);

    vector<CFeature> featin;
    vector<string> names;

    matf data;
    CFeature::OpenFromFileBlockwise("bw.txt",featin,names,data,(15*15));

//    for(size_t i=0; i<names.size(); i++) {
//        cout << featin[i] << endl;
//        cout << names[i] << endl;

//    }

    //cout << data << endl;

    return 0;


}

