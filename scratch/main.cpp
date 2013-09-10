#include <iostream>
#include <fstream>
#include "darray.h"
#include "descriptor.h"


#include "feature.h"
#include "iter.h"
#include "splinecurve.h"
#include "factor.h"
#include "types.h"
#include "trafo.h"
#include "cam.h"
#include "image.h"

#include <omp.h>
#include <stdio.h>
#include <algorithm>

#include <QImage>

using namespace std;
using namespace R4R;
using namespace cv;

/*! \brief Reads and writes feature point descriptors to file.
*/
bool rw_descriptors();

/*! \brief Finds a minimum-norm solution of over-determined linear
 * system of equations.
*/
bool solve_linear_system();

/*! \brief Loads an image and compute its gradient.
 * \param[in] filename image to load
*/
bool gradients(const char* filename);

/*! \brief Computes parallel projections into image plane of
 *  a pinhole camera.
 * \param[in] filename file containing the camera parameters
*/
bool projection(const char* filename);

int main()
{



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


bool solve_linear_system() {

    mat M(0,0);
    CPreconditioner<mat,vec,double> precond = CPreconditioner<mat,vec,double>(M);
    CIterativeSolver<mat,vec,double> solver = CIterativeSolver<mat,vec,double>(precond,50,1e-10,false);

    mat A(10,10);
    vec b(10);
    vec x(10);
    A.Rand(0,1);
    b.Rand(0,1);

    solver.CGLS(A,b,x);
    cout << x << endl;

    return 0;

}

bool gradients(const char* filename) {

    QImage qimg;
    qimg.load(QString(filename));

    CRGBImage img(qimg);

    mat Iu(img.NRows(),img.NCols()), Iv(img.NRows(),img.NCols());

    for(size_t i=0; i<Iu.NRows(); i++) {

        for(size_t j=0; j<Iv.NCols(); j++) {

            vec2 p = {(double)j,(double)i};
            vec3 gradu = img.Gradient(p,0);
            vec3 gradv = img.Gradient(p,1);

            Iu(i,j) = gradu.Get(0);
            Iv(i,j) = gradv.Get(0);


        }

    }

    Iu.WriteToFile("Iu.txt");
    Iv.WriteToFile("Iv.txt");

    return 0;

}

bool projection(const char* filename) {

    CRigidMotion<double,3> F;
    CPinholeCam cam;

    CView<double> myview(cam,F);
    myview.OpenFromFile(filename);
    cout << myview << endl;

    const CRigidMotion<double,3>& Fr = myview.GetTransformation();
    mat M = mat(Fr);
    cout << M << endl;

    CVector<size_t,2> sizes = cam.GetSize();
    cout << sizes.Get(0) << " " << sizes.Get(1) << endl;
    mat K = cam.GetProjectionMatrix();
    cout << K << endl;

    vec3 x={1,4,2};
    vec2 y = myview.Project(x);
    cout << y << endl;

    vector<vec3> in;
    for(size_t i=0; i<1000000; i++)
        in.push_back(x);

    double t0 = omp_get_wtime();
    vector<vec2> out = myview.Project(in);
    double t1 = omp_get_wtime();
    cout << "Wall clock time elapsed: " << t1-t0 << " s" << endl;

    cout << out[0] << endl;

    return 0;

}
