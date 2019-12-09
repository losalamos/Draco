//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/draco_info.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "draco_info.hh"
#include "c4/config.h"
#include "diagnostics/config.h"
#include "ds++/DracoStrings.hh"
#include "ds++/Release.hh"
#include "ds++/UnitTest.hh"
#include <algorithm> // tolower
#include <iostream>
#include <iterator>
#include <sstream>

namespace rtt_diagnostics {

//---------------------------------------------------------------------------//
//! Constructor
DracoInfo::DracoInfo(void)
    : release(rtt_dsxx::release()), copyright(rtt_dsxx::copyright()),
      contact("For information, send e-mail to draco@lanl.gov."),
      build_type(rtt_dsxx::string_toupper(CMAKE_BUILD_TYPE)),
      library_type("static"), system_type("Unknown"), site_name("Unknown"),
      cuda(false), mpi(false), mpirun_cmd(""), openmp(false),
      diagnostics_level("disabled"), diagnostics_timing(false),
      cxx(CMAKE_CXX_COMPILER), cxx_flags(CMAKE_CXX_FLAGS), cc(CMAKE_C_COMPILER),
      cc_flags(CMAKE_C_FLAGS), fc("none"), fc_flags("none") {
#ifdef DRACO_SHARED_LIBS
  library_type = "Shared";
#endif
#ifdef CMAKE_SYSTEM_NAME
  system_type = CMAKE_SYSTEM_NAME_STRING;
#endif
#ifdef SITENAME
  site_name = SITENAME;
#endif
#ifdef HAVE_CUDA
  cuda = true;
#endif
#ifdef C4_MPI
  mpi = true;
  mpirun_cmd = std::string(MPIEXEC_EXECUTABLE) + std::string(" ") +
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
  if (build_type == std::string("Release")) {
    cxx_flags += CMAKE_CXX_FLAGS_RELEASE;
    cc_flags += CMAKE_C_FLAGS_RELEASE;
  } else if (build_type == std::string("RelWithDebInfo")) {
    cxx_flags += CMAKE_CXX_FLAGS_RELWITHDEBINFO;
    cc_flags += CMAKE_C_FLAGS_RELWITHDEBINFO;
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
void print_text_with_word_wrap(std::string const &longstring,
                               size_t const indent_column,
                               size_t const max_width, std::ostringstream &msg,
                               std::string const &delimiters = " ") {
  std::vector<std::string> const tokens =
      rtt_dsxx::tokenize(longstring, delimiters);
  std::string const delimiter(delimiters.substr(0, 1));
  size_t i(indent_column);
  for (auto item : tokens) {
    if (i + item.length() + 1 > max_width) {
      msg << "\n" << std::string(indent_column, ' ');
      i = indent_column;
    }
    msg << item;
    if (item != tokens.back())
      msg << delimiter;
    i += item.length() + 1;
  }
}

//---------------------------------------------------------------------------//
std::string DracoInfo::fullReport(void) const {
  using std::cout;
  using std::endl;

  std::ostringstream infoMessage;

  // Create a list of features for DbC
  std::vector<std::string> dbc_info;
  dbc_info.push_back("Insist");
#ifdef REQUIRE_ON
  dbc_info.push_back("Require");
#endif
#ifdef CHECK_ON
  dbc_info.push_back("Check");
#endif
#ifdef ENSURE_ON
  dbc_info.push_back("Ensure");
#endif
#if DBC & 8
  dbc_info.push_back("no-throw version");
#endif
#if DBC & 16
  dbc_info.push_back("check-deferred version");
#endif

  // Print version and copyright information to the screen:
  infoMessage << briefReport();

  // Build Information
  //------------------

  infoMessage << "Build information:"
              << "\n    Build type        : " << build_type
              << "\n    Library type      : " << library_type
              << "\n    System type       : " << system_type
              << "\n    Site name         : " << site_name
              << "\n    CUDA support      : " << (cuda ? "enabled" : "disabled")
              << "\n    MPI support       : "
              << (mpi ? "enabled" : "disabled (c4 scalar mode)");

  if (mpi)
    infoMessage << "\n      mpirun cmd      : " << mpirun_cmd;

  infoMessage << "\n    OpenMP support    : "
              << (openmp ? "enabled" : "disabled")
              << "\n    Design-by-Contract: " << DBC << ", features = ";
  std::copy(dbc_info.begin(), dbc_info.end() - 1,
            std::ostream_iterator<std::string>(infoMessage, ", "));
  infoMessage << dbc_info.back();
  infoMessage << "\n    Diagnostics       : " << diagnostics_level
              << "\n    Diagnostics Timing: "
              << (diagnostics_timing ? "enabled" : "disabled");

  // Compilers and Flags
  size_t const max_width(80);
  size_t const hanging_indent(std::string("    CXX Compiler      : ").length());
  infoMessage << "\n    CXX Compiler      : ";
  print_text_with_word_wrap(cxx, hanging_indent, max_width, infoMessage, "/");
  infoMessage << "\n    CXX_FLAGS         : ";
  print_text_with_word_wrap(cxx_flags, hanging_indent, max_width, infoMessage);
  infoMessage << "\n    C Compiler        : ";
  print_text_with_word_wrap(cc, hanging_indent, max_width, infoMessage, "/");
  infoMessage << "\n    C_FLAGS           : ";
  print_text_with_word_wrap(cc_flags, hanging_indent, max_width, infoMessage);
  infoMessage << "\n    Fortran Compiler  : ";
  print_text_with_word_wrap(fc, hanging_indent, max_width, infoMessage, "/");
  infoMessage << "\n    Fortran_FLAGS     : ";
  print_text_with_word_wrap(fc_flags, hanging_indent, max_width, infoMessage);

  infoMessage << "\n" << endl;

  return infoMessage.str();
}

//---------------------------------------------------------------------------//
std::string DracoInfo::briefReport(void) const {
  std::ostringstream infoMessage;

  // Print version and copyright information to the screen:
  infoMessage << "\n";
  print_text_with_word_wrap(release, 5, 80, infoMessage, ";");
  infoMessage << "\n\n" << copyright << "\n" << contact << "\n" << std::endl;
  return infoMessage.str();
}

//---------------------------------------------------------------------------//
//! extract the single-line version info from release and return it
std::string DracoInfo::versionReport(void) const {
  std::ostringstream infoMessage;
  print_text_with_word_wrap(release, 5, 80, infoMessage, ";");
  infoMessage << "\n" << std::endl;
  return infoMessage.str();
}

} // end namespace rtt_diagnostics

//---------------------------------------------------------------------------//
// end of draco_info.cc
//---------------------------------------------------------------------------//
