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

#include <OpenEXR/ImfArray.h>
#include <OpenEXR/ImfRgbaFile.h>
#include <mex.h>

using namespace Imf;
using namespace Imath;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // check arguments
    if(nrhs != 1 || mxIsChar(prhs[0])!=1)
        mexErrMsgTxt("Wrong number of input arguments.");

    if(nlhs != 4)
        mexErrMsgTxt("Wrong number of output arguments.");

    // filename 
    char* fname = mxArrayToString(prhs[0]);
	
    try {
		
        RgbaInputFile file(fname);

        Box2i dw = file.dataWindow();

        size_t w = dw.max.x - dw.min.x + 1;
        size_t h = dw.max.y - dw.min.y + 1;

        Array2D<Rgba> in(h,w);

        file.setFrameBuffer(&in[0][0]-dw.min.x-dw.min.y*w,1,w);
        file.readPixels(dw.min.y, dw.max.y);
       
        // convert dimensions
		int dims[2];
		dims[0] = h;
		dims[1] = w;

        // allocate memory in matlab workspace
        plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        plhs[3] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
        
        // get ptrs to memory
        double* r = mxGetPr(plhs[0]);
        double* g = mxGetPr(plhs[1]);
        double* b = mxGetPr(plhs[2]);
        double* a = mxGetPr(plhs[3]);

		for (size_t i=0; i<h; i++) {
            
			for (size_t j=0; j<w; j++) {	
                
                size_t k = j*h + i;

                r[k] = (float)in[i][j].r;
                g[k] = (float)in[i][j].g;
                b[k] = (float)in[i][j].b;
                a[k] = (float)in[i][j].a;
                
			}
            
		}
        
	}
    catch (const exception& e) {
		
        mxFree(fname);
		mexErrMsgTxt(e.what());
	
    }

}
