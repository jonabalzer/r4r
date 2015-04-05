### Features


* Linear algebra
    * Sparse and dense matrix data structures
    * Iterative linear solvers
        * CG
        * CGLS
    * Preconditioners
       * Jacobi
       * SSOR
* LAPACK wrappers
* Optimization
    * Levenberg-Marquardt
    * Split-Bregman
    * Reweighted least-squares
    * PEGASOS SVM solver
    * Kalman filter
* Tracking
    * Feature point data structure
    * Track administration
    * TST
* Geometry
    * Camera models
    * Coordinate transformations
    * B-splines in arbitrary dimensions
* Image descriptors
    * BRIEF
    * HOG
    * Raw image (with different normalizations)
    * Import/export
* Example applications
    * Robust estimation
    * CLAM (real-time structure from motion)
    * Multiview descriptor aggregation
    * TV image denoising
* Python bindings
    * File I/O
    * B-splines
* Documentation
* Unit testing

### Quickstart

#### Dependencies

- [OpenCV 3.0 (!)](http://opencv.org/)
- [LAPACK](http://www.netlib.org/lapack/)
- [QT](http://qt-project.org/) (only required by examples with GUI)
- [OpenMesh](www.openmesh.org/‎) (required by parts of the reconstruction library)
- [FFTW](http://www.fftw.org/) (optional)

An automatic build requires the following items:

- [qmake](http://qt-project.org/)
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (documentation only)
- [Latex](http://www.tug.org/texlive/) (documentation only)
- [GraphViz](http://www.graphviz.org/) (documentation only)

#### Build instructions (LINUX)


Make sure that aforementioned dependencies are installed in default locations. __In particular, it is highly recommended to build the latest OpenCV version (i.e., the trunk) from source.__ To build R4R, clone this repository by typing
```
git clone https://github.com/jonabalzer/r4r.git
```
This will automatically create a subfolder `r4r`. Change into this directory
```
cd r4r
```
and create a directory for the out-of-core build, say
```
mkdir build
```
Call qmake from the build directory:
```
cd build
qmake ..
```
If you want to include the example applications, you need to set the `HAVE_EXAMPLES` variable:
```
qmake .. "HAVE_EXAMPLES=1"
```
Start the build process by
```
make
```
The documentation is created via
```
make doc
```
If you plan to use the core libraries outside of Qt Creator, run
```
sudo make install
```
This will install libraries and header files into the appropriate system paths (which probably requires root privileges).
