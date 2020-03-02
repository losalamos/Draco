Draco
----------------

[![Linux Build Status](https://travis-ci.org/lanl/Draco.svg?branch=develop)](https://travis-ci.org/lanl/Draco)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/yp8r9jxl2gc9n1fs/branch/develop?svg=true)](https://ci.appveyor.com/project/lanl/Draco)
[![codecov.io](https://codecov.io/github/lanl/Draco/coverage.svg?branch=develop)](https://codecov.io/github/lanl/Draco/branch/develop)
[![Latest Version](https://img.shields.io/github/release/lanl/draco.svg?style=flat-square)](https://github.com/lanl/Draco/releases)
[![PyPI](https://img.shields.io/pypi/l/Django.svg)](https://github.com/lanl/Draco/blob/develop/LICENSE.md)

Draco is an object-oriented component library geared towards numerically
intensive, radiation (particle) transport applications built for parallel
computing hardware.  It consists of semi-independent packages and a robust build
system.  The packages in draco provide a set of components that can be used by
multiple clients to build transport codes.  The build system can also be
extracted for use in clients.

To clone draco:

    $ git clone https://github.com/lanl/Draco.git

To get started, please see [Development - Quick Start Guide](https://github.com/lanl/Draco/wiki/Development---Quick-Start)
in the wiki. Pull requestes must satisfy the requirements listed in
the [Style Guide](https://github.com/lanl/Draco/wiki/Style-Guide).

Synopsis of Active Draco Packages
---------------------------------

* c4 - A communications library for message passing interfaces (MPI).
  * For builds without MPI, all communication commnds will be no-op functions
    (`DRACO_C4={MPI;SCALAR}`).
* cdi - Access to material data. The Common Data Interface (CDI) specifies a
  common abstraction for objects and libraries that return material data
  (opacities, atomic cross sections, equation-of-state data, etc.)
* cdi_analytic - Analytic models for physical data
* cdi_eospac - Equation-of-State data
  * These classes wrap the EOSPAC6 libraries that read sesame files; Commonly
    used to access gray opacity data and heat capacities.
   (`EOSPAC_LIB_DIR=<path>`).
* cdi_ipcress - Gray and multigroup opacities
  * The classes in this component will read and parse opacity values from an
    IPCRESS file produced by TOPS.
* compton - Provides access to Compton scattering models and data as provided
  by the CSK library.
* device - Wrapper for heterogeneous device communication
  * The classes in this component provide access to DaCS (deprecated) and CUDA
    calls for use on heterogeneous architecture platforms (GPU machines).
* diagnostics - CPP macros that are activated at compile time that can provide
  additional diagnostic verbosity during calculations.
* ds++ - Basic services and data structures library.
  * Array containers, assertion and Design-by-Contract, access to low level OS
    functions, file and path manipulation, unit test system, etc.
* fit - Least squares fitting routines.
* fpe_trap - Catch IEEE floating point exceptions
* FortranCheck - Test Fortran compatibility and interoperability
  * The examples in this component will demonstrate if the Fortran compiler is working; if Fortran/C interlanguage linking/running is working and sample ISO_C_BINDING calls.
* lapack_wrap - C++ wrapper for BLAS and LAPACK.
* linear - direct solvers for small linear systems of equations.
* mesh - Encapsulate mesh definition and accessors.
* mesh_element - defines fundamental mesh element types used by mesh,
  meshReaders and RTT_Format_Reader.
* meshReaders - Read RTT format mesh files.
* min - Optimization routines. Find the minimum of a function.
* norms - Calculate norms for data sets.
* ode - Ordinary differential equation solvers (e.g.: Runge-Kutta).
* parser - Lexical file/input parser.
* plot2d - C++ interface to XMGrace 2-dimensional plotting.
* quadrature - access to angular integration functions and related data.
* rng - A random number generator component
  * The primary set of functions provided by this component were derived from
    Random123 (https://www.deshawresearch.com/downloads/download_random123.cgi)
    random number library.  A few additional random number generators are also
    provided.
* roots - Root finding algorithms
* RTT_Format_Reader - meshReaders implementation for RTT format mesh files or
  input-streams.
* shared_lib - Dynamically load/unload shared object libraries via dl load
* special_functions - Specialized math functions like factorial and Dirac delta.
* traits - A traits class used by viz.
* units - Provides encapsulated unit systems, functions to convert between unit
  systems and physical constants.
* VendorChecks - A testing component to ensure that discovered and used vendor
  libraries behave as expected.
* viz - Generates Ensight files for data visualization.

Authors
----------------
Many thanks go to Draco's [contributors](https://github.com/lanl/Draco/graphs/contributors).

Draco was originally written by staff from Los Alamos's [CCS-2 Computational Physics and Methods Group](http://www.lanl.gov/org/padste/adtsc/computer-computational-statistical-sciences/computational-physics-methods/index.php):

> *CCS-2 Draco Team:* Kelly G. Thompson, Kent G. Budge, Ryan T. Wollaeger,
>   James S. Warsa, Alex R. Long, Kendra P. Long, Jae H. Chang,
>   Matt A. Cleveland, Andrew T. Till, and Tim Kelley.

> *Prior Contributers:* Gabriel M. Rockefeller, Allan B. Wollaber,
>   Lori A. Pritchett-Sheats, Rob B. Lowrie, Paul W. Talbot, Katherine J. Wang,
>   Peter Ahrens, Daniel Holladay, Jeff D. Densmore, Massimiliano Rosa,
>   Todd J. Urbatsch, Jeff Furnish, John McGhee, Kris C. Garrett, Mike Buksas,
>   Nick Myers, Paul Henning, Randy Roberts, Seth Johnson, Todd Adams, and
>   Tom Evans.

Release
----------------

Draco is released under the BSD 3-Clause License. For more details see the
[LICENSE file](https://github.com/lanl/Draco/blob/develop/LICENSE.md).
