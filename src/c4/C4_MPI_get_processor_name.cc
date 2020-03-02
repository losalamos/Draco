//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/C4_MPI_get_processor_name.cc
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI implementation.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//
// $Id: C4_MPI_get_processor_name.cc 7388 2015-01-22 16:02:07Z kellyt $
//----------------------------------------------------------------------------//

#include "c4/config.h"
#include <string>

#ifdef C4_MPI

#include "C4_Functions.hh"

namespace rtt_c4 {
//----------------------------------------------------------------------------//
// get_processor_name
//----------------------------------------------------------------------------//
std::string get_processor_name() {
  int namelen(0);
  char processor_name[DRACO_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(processor_name, &namelen);
  std::string pname(processor_name);
  Ensure(pname.size() == static_cast<size_t>(namelen));
  return pname;
}

} // end namespace rtt_c4

#endif // C4_MPI

//----------------------------------------------------------------------------//
// end of C4_MPI_get_processor_name.cc
//----------------------------------------------------------------------------//
