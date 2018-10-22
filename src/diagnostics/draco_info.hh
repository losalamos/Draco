//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/draco_info.hh
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_diagnostics_draco_info_hh
#define rtt_diagnostics_draco_info_hh

#include "ds++/config.h"
#include <string>
#include <vector>

namespace rtt_diagnostics {

//===========================================================================//
/*!
 * \class DracoInfo
 * \brief Store and present basic information about the current draco build.
 *
 * The constructed string will take this form:
 * \verbatim
 * Draco-6_5_20121113, build date 2012/11/13; build type: DEBUG; DBC: 7
 *
 * Draco Contributers:
 *     Kelly G. Thompson, Kent G. Budge, Tom M. Evans,
 *     Rob Lowrie, B. Todd Adams, Mike W. Buksas,
 *     James S. Warsa, John McGhee, Gabriel M. Rockefeller,
 *     Paul J. Henning, Randy M. Roberts, Seth R. Johnson,
 *     Allan B. Wollaber, Peter Ahrens, Jeff Furnish,
 *     Paul W. Talbot, Jae H. Chang, and Benjamin K. Bergen.
 *
 * Copyright (C) 2016 LANS, LLC
 * Build information:
 *     Library type   : shared
 *     System type    : Linux
 *     CUDA support   : disabled
 *     MPI support    : enabled
 *       mpirun cmd   : mpirun --mca mpi_paffinity_alone 0 -np
 *     OpenMPI support: enabled
 *     Diagnostics    : disabled
 *     Diagnostics Timing: disabled
 * \endverbatim
 */
//===========================================================================//
class DLL_PUBLIC_diagnostics DracoInfo {
public:
  // IMPLELEMENTATION
  // ================

  // Constructors
  // ------------
  DracoInfo(void);

  // Actions
  // -------

  /*! \brief Construct an information message that includes Draco's version,
   * copyright and basic build parameters. */
  std::string fullReport(void);

  //! Version and Copyright only
  std::string briefReport(void);

  //! Version only
  std::string versionReport(void);

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

//---------------------------------------------------------------------------//
// end of draco_info.hh
//---------------------------------------------------------------------------//
