//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   diagnostics/draco_info.hh
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef rtt_diagnostics_draco_info_hh
#define rtt_diagnostics_draco_info_hh

#include "ds++/config.h"
#include <string>
#include <vector>

namespace rtt_diagnostics {

//============================================================================//
/*!
 * \class DracoInfo
 * \brief Store and present basic information about the current draco build.
 *
 * The constructed string will take this form:
 * \verbatim
Draco-6_25_20181114, build date 2018/11/14;build type: Debug;DBC: 7;
     DRACO_DIAGNOSTICS: 0

CCS-2 Draco Team: Kelly G. Thompson, Kent G. Budge, Ryan T. Wollaeger,
     James S. Warsa, Alex R. Long, Kendra P. Keady, Jae H. Chang,
     Matt A. Cleveland, Andrew T. Till, Tim Kelley, and Kris C. Garrett.
Prior Contributers: Jeff D. Densmore, Gabriel M. Rockefeller,
     Allan B. Wollaber, Rob B. Lowrie, Lori A. Pritchett-Sheats,
     Paul W. Talbot, and Katherine J. Wang.

Copyright (C) 2016-2020 Triad National Security, LLC. (LA-CC-16-016)

For information, send e-mail to draco@lanl.gov.

Build information:
    Build type        : DEBUG
    Library type      : Shared
    System type       : Linux
    Site name         : ccscs3
    CUDA support      : disabled
    MPI support       : enabled
      mpirun cmd      : /scratch/.../bin/mpiexec -n <N> -bind-to none
    OpenMP support    : enabled
    Design-by-Contract: 7, features = Insist, Require, Check, Ensure
    Diagnostics       : 0
    Diagnostics Timing: disabled
    CXX Compiler      : scratch/vendors/spack.20180425/opt/spack/
                        linux-rhel7-x86_64/gcc-4.8.5/
                        gcc-8.1.0-3c5hjkqndywdp3w2l5vts62xlllrsbtq/bin/g++
    CXX_FLAGS         : -Wcast-align -Wpointer-arith -Wall -pedantic
                        -Wno-expansion-to-defined -Wnarrowing -march=native
                        -fopenmp -Werror
    C Compiler        : scratch/vendors/spack.20180425/opt/spack/
                        linux-rhel7-x86_64/gcc-4.8.5/
                        gcc-8.1.0-3c5hjkqndywdp3w2l5vts62xlllrsbtq/bin/gcc
    C_FLAGS           : -Wcast-align -Wpointer-arith -Wall -pedantic
                        -Wno-expansion-to-defined -Wnarrowing -march=native
                        -fopenmp -Werror
    Fortran Compiler  : scratch/vendors/spack.20180425/opt/spack/
                        linux-rhel7-x86_64/gcc-4.8.5/
                        gcc-8.1.0-3c5hjkqndywdp3w2l5vts62xlllrsbtq/bin/gfortran
    Fortran_FLAGS     : -ffree-line-length-none -cpp -march=native -fopenmp

 * \endverbatim
 */
//============================================================================//
class DracoInfo {
public:
  // IMPLELEMENTATION
  // ================

  // Constructors
  // ------------
  DracoInfo(void);

  // Actions
  // -------

  /*! \brief Construct an information message that includes Draco's version,
   *         copyright and basic build parameters. */
  std::string fullReport(void) const;

  //! Version and Copyright only
  std::string briefReport(void) const;

  //! Version only
  std::string versionReport(void) const;

private:
  // DATA
  // ----

  std::string const release;
  std::string const copyright;
  std::string const contact;
  std::string const build_type;
  std::string library_type;
  std::string system_type;
  std::string site_name;
  bool cuda;
  bool mpi;
  std::string mpirun_cmd;
  bool openmp;
  std::string diagnostics_level;
  bool diagnostics_timing;
  std::string cxx;
  std::string cxx_flags;
  std::string cc;
  std::string cc_flags;
  std::string fc;
  std::string fc_flags;
};

} // end namespace rtt_diagnostics

#endif //  rtt_diagnostics_draco_info_hh

//----------------------------------------------------------------------------//
// end of draco_info.hh
//----------------------------------------------------------------------------//
