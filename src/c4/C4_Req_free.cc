//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Req_free.cc
 * \author Thomas M. Evans, Geoffrey Furnish
 * \date   Thu Jun  2 09:54:02 2005
 * \brief  C4_Req member definitions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: C4_Req_free.cc 7388 2015-01-22 16:02:07Z kellyt $
//---------------------------------------------------------------------------//

#include "C4_Req.hh"

namespace rtt_c4 {

void C4_ReqRefRep::free() {
#ifdef C4_MPI
  if (assigned) {
    MPI_Cancel(&r);
    MPI_Request_free(&r);
  }
#endif
  clear();
}

//---------------------------------------------------------------------------//
//! Return the number of items returned on the last complete operation.
//---------------------------------------------------------------------------//

unsigned C4_ReqRefRep::count() {
#ifdef C4_MPI
  int count;
  MPI_Get_count(&s, MPI_CHAR, &count);
  return count;
#else
  return 0;
#endif
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of C4_Req_free.cc
//---------------------------------------------------------------------------//
