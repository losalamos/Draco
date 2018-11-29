Draco
----------------

[![Linux Build Status](https://travis-ci.org/lanl/Draco.svg?branch=develop)](https://travis-ci.org/lanl/Draco)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/yp8r9jxl2gc9n1fs/branch/develop?svg=true)](https://ci.appveyor.com/project/lanl/Draco)
[![codecov.io](https://codecov.io/github/lanl/Draco/coverage.svg?branch=develop)](https://codecov.io/github/lanl/Draco/branch/develop)
[![Latest Version](https://img.shields.io/github/release/lanl/draco.svg?style=flat-square)](https://github.com/lanl/Draco/releases)
[![PyPI](https://img.shields.io/pypi/l/Django.svg)](https://github.com/lanl/Draco/blob/develop/LICENSE.md)

Draco is an object-oriented component library geared towards
numerically intensive, radiation (particle) transport applications
built for parallel computing hardware.  It consists of
semi-independent packages and a robust build system.  The packages in
draco provide a set of components that can be used by multiple clients
to build transport codes.  The build system can also be extracted for
use in clients.

To clone draco:

    $ git clone https://github.com/lanl/Draco.git

To get started, please see [Development - Quick Start Guide](https://github.com/lanl/Draco/wiki/Development---Quick-Start)
in the wiki. Pull requestes must satisfy the requirements listed in
the [Style Guide](https://github.com/lanl/Draco/wiki/Style-Guide).

Synopsis of Active Draco Packages
---------------------------------

* RTT_Format_Reader - mesh container and services
* c4 - Communications library that wraps MPI
* cdi - Get access to material data
  * The Common Data Interface specifies an common abstraction for objects and libraries that return material data (opacities, atomic cross sections, equation-of-state data, etc.)
* cdi_analytic - Analytic models for physical data
* cdi_eospac - Equation-of-State data
  * This class wraps the EOSPAC6 libraries. The interface inherits from CDI/EOS.hh and can be used as a stand alone component or as a plug-in to CDI.  To make the component available you must set in your environment or in the CMakeCache.txt the variable EOSPAC_LIB_DIR to the path to libeospac.a
* cdi_ipcress - Gray and multigroup opacities
  *  The classes in this component will read and parse opacity values from an IPCRESS file produced by TOPS.  The component presents this data through various accessors.  The interface inherits from CDI/Opacity.hh and can be used as a stand alone component or as a plug-in to CDI.
* device - Wrapper for heterogeneous device communication
  * The classes in this component provide access to DaCS and CUDA calls for use on heterogeneous architecture platforms (Roadrunner or GPU machines).
* diagnostics - CPP macros that are activated at compile time that can provide additional diagnostic verbosity during calculations.
* ds++ - Basic services
  * Array containers, assertion and Design-by-Contract, Ffle system information, path manipulation, unit test system, etc.
* fit - Least squares fit
* fpe_trap - Catch IEEE floating point exceptions
* FortranCheck - Test Fortran compatibility and interoperability
  * The examples in this component will demonstrate if the Fortran compiler is working; if Fortran/C interlanguage linking/running is working and sample ISO_C_BINDING calls.
* lapack_wrap - BLAS and LAPACK functionality for C++
* linear - linear equations
  * Routines to solve small systems of linear equations.
* meshReaders - Read RTT format mesh files.
* mesh - Encapsulate mesh definition and accessors.
* mesh_element - Provide a description of unstructured mesh. Used by meshReaders and RTT_Format_Reader.
* min - Optimizaiton routines. Find the minimum of a function.
* norms - Calculate norms for data sets.
* ode - Ordinary differential equation solvers.
* parser - Lexical file parser.
* plot2d - Generate GNU plot 2-dimensional plots.
* quadrature - Get access to quadrature data
  * Provides a creator class to generate quadrature objects.  Quadrature objects encapsulate discrete ordinates angular differencing methods and data. Also provide service functions that are commonly used.
* rng - A random number generator component
  * The primary set of functions provided by this component were derived from Random123 (https://www.deshawresearch.com/downloads/download_random123.cgi) random number library.  A few additional random number generators are also provided.
* roots - Root solvers
* shared_lib - Dynamically load/unload shared object libraries via dl load
* special_functions - Specialized math functions like factorial and Dirac delta.
* timestep - An object-oriented class that encapsulates a time step controller.
* traits - A traits class used by viz.
* units:
  * Provides units standardization for Draco. Contains Units class for representing arbitrary systems and converting quantities to SI units.
* viz - Generates Ensight files for data vizualization.

Authors
----------------
Many thanks go to Draco's [contributors](https://github.com/lanl/Draco/graphs/contributors).

Draco was originally written by staff from Los Alamos's [CCS-2 Computational Physics and Methods Group](http://www.lanl.gov/org/padste/adtsc/computer-computational-statistical-sciences/computational-physics-methods/index.php):

> *CCS-2 Draco Team:* Kelly G. Thompson, Kent G. Budge, Ryan T. Wollaeger,
> James S. Warsa, Alex R. Long, Kendra P. Keady, Jae H. Chang,
> Matt A. Cleveland, Andrew T. Till, Tim Kelley, and Kris C. Garrett.

> *Prior Contributers:* Jeff D. Densmore, Gabriel M. Rockefeller,
> Allan B. Wollaber, Rob B. Lowrie, Lori A. Pritchett-Sheats,
> Paul W. Talbot, and Katherine J. Wang.

Release
----------------

Draco is released under the BSD 3-Clause License. For more details see the
[LICENSE file](https://github.com/lanl/Draco/blob/develop/LICENSE.md).

LANL code designation: `LA-CC-16-016`
