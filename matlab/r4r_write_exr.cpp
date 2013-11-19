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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // check arguments
    if(nrhs != 5 || mxIsChar(prhs[0])!=1)
        mexErrMsgTxt("Wrong number of input arguments.");

    if(nlhs != 0)
        mexErrMsgTxt("Wrong number of output arguments.");

    // filename 
    char* fname = mxArrayToString(prhs[0]);

    // get ptrs to memory
    double* r = mxGetPr(prhs[1]);
    double* g = mxGetPr(prhs[2]);
    double* b = mxGetPr(prhs[3]);
    double* a = mxGetPr(prhs[4]);
    
    // get size
    int height = mxGetM(prhs[1]);
    int width = mxGetN(prhs[1]);
	
    // create openexr array
    Array2D<Rgba> out(height,width);
    
    try {
		      
		for (size_t i=0; i<height; i++) {
            
			for (size_t j=0; j<width; j++) {	
                
                size_t k = j*height + i;
                
                Rgba val;
                val.r = half(r[k]);
                val.g = half(g[k]);
                val.b = half(b[k]);
                val.a = half(a[k]);
                                
                out[i][j] = val;
                
			}
            
		}
        
         RgbaOutputFile file(fname, width, height, WRITE_RGBA);
         file.setFrameBuffer (&out[0][0],1,width);
         file.writePixels (height);
        
	}
    catch(const exception& e) {
		
        mxFree(fname);
		mexErrMsgTxt(e.what());
	
    }

}
