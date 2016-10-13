//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/draco_info.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstScalarUnitTest.cc 6864 2012-11-08 01:34:45Z kellyt $
//---------------------------------------------------------------------------//

#include "draco_info.hh"
#include "c4/config.h"
#include "diagnostics/config.h"
#include "ds++/Release.hh"
#include "ds++/UnitTest.hh"
#include <algorithm> // tolower
#include <iostream>
#include <sstream>

namespace rtt_diagnostics {

//---------------------------------------------------------------------------//
//! Constructor
DracoInfo::DracoInfo(void)
    : release(rtt_dsxx::release()), copyright(rtt_dsxx::copyright()),
      contact("For information, send e-mail to draco@lanl.gov."),
      build_type(normalizeCapitalization(CMAKE_BUILD_TYPE)),
      library_type("static"), system_type("Unknown"), site_name("Unknown"),
      cuda(false), mpi(false), mpirun_cmd(""), openmp(false),
      diagnostics_level("disabled"), diagnostics_timing(false), cxx11(false),
      cxx11_features(), cxx(CMAKE_CXX_COMPILER), cxx_flags(CMAKE_CXX_FLAGS),
      cc(CMAKE_C_COMPILER), cc_flags(CMAKE_C_FLAGS), fc("none"),
      fc_flags("none") {
#ifdef DRACO_SHARED_LIBS
  library_type = "Shared";
#endif
#if DRACO_UNAME == Linux
  system_type = "Linux";
#endif
#ifdef SITENAME
  site_name = SITENAME;
#endif
#ifdef HAVE_CUDA
  cuda = true;
#endif
#ifdef C4_MPI
  mpi = true;
  mpirun_cmd = std::string(MPIEXEC) + std::string(" ") +
               std::string(MPIEXEC_NUMPROC_FLAG) + std::string(" <N> ");
#ifdef MPIEXEC_POSTFLAGS
  mpirun_cmd += std::string(MPIEXEC_POSTFLAGS);
#endif
#endif
#ifdef OPENMP_FOUND
  openmp = true;
#endif
#ifdef DRACO_DIAGNOSTICS
  std::ostringstream msg;
  msg << DRACO_DIAGNOSTICS;
  diagnostics_level = msg.str();
#endif
#ifdef DRACO_TIMING
  diagnostics_timing = true;
#endif
  cxx11 = true;
  cxx11_features = rtt_dsxx::UnitTest::tokenize(CXX11_FEATURE_LIST, ";", false);
  if (build_type == std::string("Release")) {
    cxx_flags += CMAKE_CXX_FLAGS_RELEASE;
    cc_flags += CMAKE_C_FLAGS_RELEASE;
  } else if (build_type == std::string("Debug")) {
    cxx_flags += CMAKE_CXX_FLAGS_DEBUG;
    cc_flags += CMAKE_C_FLAGS_DEBUG;
  }
#ifdef CMAKE_Fortran_COMPILER
  fc = CMAKE_Fortran_COMPILER;
  fc_flags = CMAKE_Fortran_FLAGS;
  if (build_type == std::string("Release"))
    fc_flags += CMAKE_Fortran_FLAGS_RELEASE;
  else if (build_type == std::string("Debug"))
    fc_flags += CMAKE_Fortran_FLAGS_DEBUG;
#endif
}

//---------------------------------------------------------------------------//
std::string DracoInfo::fullReport(void) {
  using std::cout;
  using std::endl;

  std::ostringstream infoMessage;

  // Print version and copyright information to the screen:
  infoMessage << briefReport();

  // Build Information
  //------------------

  infoMessage << "Build information:"
              << "\n    Build type     : " << build_type
              << "\n    Library type   : " << library_type
              << "\n    System type    : " << system_type
              << "\n    Site name      : " << site_name
              << "\n    CUDA support   : " << (cuda ? "enabled" : "disabled")
              << "\n    MPI support    : "
              << (mpi ? "enabled" : "disabled (c4 scalar mode)");

  if (mpi)
    infoMessage << "\n      mpirun cmd   : " << mpirun_cmd;

  infoMessage << "\n    OpenMP support : " << (openmp ? "enabled" : "disabled")
              << "\n    Diagnostics    : " << diagnostics_level
              << "\n    Diagnostics Timing: "
              << (diagnostics_timing ? "enabled" : "disabled");

  // C++11

  infoMessage << "\n    C++11 Support  : " << (cxx11 ? "enabled" : "disabled");

  if (cxx11) {
    infoMessage << "\n      Feature list : ";
    for (size_t i = 0; i < cxx11_features.size(); ++i)
      if (i == 0)
        infoMessage << cxx11_features[i];
      else
        infoMessage << "\n                     " << cxx11_features[i];
  }

  // Compilers and Flags

  infoMessage << "\n    CXX Compiler      : " << cxx
              << "\n    CXX_FLAGS         : " << cxx_flags
              << "\n    C Compiler        : " << cc
              << "\n    C_FLAGS           : " << cc_flags
              << "\n    Fortran Compiler  : " << fc
              << "\n    Fortran_FLAGS     : " << fc_flags;

  infoMessage << "\n" << endl;

  return infoMessage.str();
}

//---------------------------------------------------------------------------//
std::string DracoInfo::briefReport(void) {
  std::ostringstream infoMessage;

  // Print version and copyright information to the screen:
  infoMessage << "\n"
              << release << "\n\n"
              << copyright << "\n"
              << contact << "\n"
              << std::endl;

  return infoMessage.str();
}

//---------------------------------------------------------------------------//
//! extract the single-line version info from release and return it
std::string DracoInfo::versionReport(void) {
  std::ostringstream infoMessage;
  infoMessage << release << "\n" << std::endl;
  return infoMessage.str();
}

//---------------------------------------------------------------------------//
// Create a string to hold build_type and normalize the case.

std::string DracoInfo::normalizeCapitalization(std::string mystring) {
  std::transform(mystring.begin(), mystring.end(), mystring.begin(), ::tolower);
  mystring[0] = ::toupper(mystring[0]);
  return mystring;
}

} // end namespace rtt_diagnostics

//---------------------------------------------------------------------------//
// end of draco_info.cc
//---------------------------------------------------------------------------//
