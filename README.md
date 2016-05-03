Draco
----------------

Draco is an object-oriented component library geared towards
numerically intensive, radiation (particle) transport applications
built for parallel computing hardware.  It consists of
semi-independent packages and a robust build system.  The packages in
draco provide a set of components that can be used by multiple clients
to build transport codes.  The build system can also be extracted for
use in clients.

To clone draco:

    $ git clone https://github.com/losalamos/Draco.git

Compiling
----------------

Prerequisites:

* [cmake-3.5 or later](https://cmake.org/download/) 
  * cmake must be in your PATH
* A modern C++ compiler (g++, icpc, pgCC, clang++, cl)
* MPI (openMPI, mpich) 
  * mpiexec must be in your PATH
* [Random123](https://www.deshawresearch.com/downloads/download_random123.cgi)
  * export RANDOM123_INC_DIR=/some/path/Random123-1.08/include
* [Gnu Scientific Library](http://www.gnu.org/software/gsl/)
  * gsl-config must be in your PATH
* python 2X

Only needed for testing:
* [numdiff](https://www.nongnu.org/numdiff)
  * numdiff must be in your PATH

Configure:
* Use a separate build directory
```
    $ mkdir build
    $ cmake ../Draco.git
```
* For Cray PrgEnv, you need to prime the CMakeCache.txt
```
    $ cmake -C ../Draco.git/config/CracyConfig.cmake ../Draco.git
```

Optional components:

* Additional libraries will be built if certain features are available in the build environment:
  * `lapack_wrap` will be built if [LAPACK](http://www.netlib.org/lapack) or equivalent is available.
  * `cdi_eospac` will be built if [LANL's libeospac](http://www.lanl.gov/org/padste/adtsc/theoretical/physics-chemistry-materials/sesame-database.php) is available.
  * `device` will be built if the CUDA toolkit is available and the build machine has a GPU.
  * `plot2D` will be built if the Grace headers and library are found.

Build:
```
   $ make
```
Test:
```
   $ ctest
```
Install:
```
   $ make install
```

Authors
----------------
Many thanks go to Draco's [contributors](https://github.com/losalamos/Draco/graphs/contributors).

Draco was originally written by staff from Los Alamos's [CCS-2 Computational Physics and Methods Group](http://www.lanl.gov/org/padste/adtsc/computer-computational-statistical-sciences/computational-physics-methods/index.php):

> Kelly G. Thompson, Kent G. Budge, Tom M. Evans, B. Todd Adams,
> Rob Lowrie, John McGhee, Mike W. Buksas, Gabriel M. Rockefeller,
> James S. Warsa, Seth R. Johnson, Allan B. Wollaber, Randy M. Roberts,
> Jae H. Chang, Paul W. Talbot, Katherine J. Wang, Jeff Furnish,
> Matthew A. Cleveland, Benjamin K. Bergen, Paul J. Henning.

Release
----------------

Draco is released under the BSD 3-Clause License. For more details see the
LICENSE file.

LANL code designation: `LA-CC-16-016`
