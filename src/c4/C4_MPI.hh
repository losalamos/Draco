//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:16 2002
 * \brief  C4 MPI function declarations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_C4_MPI_hh
#define c4_C4_MPI_hh

#include <c4/config.h>

#ifdef C4_MPI

#include "C4_Functions.hh"
#include "ds++/Assert.hh"
#include "c4_mpi.h"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

extern MPI_Comm communicator;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

template<class Comm>
void inherit(const Comm &comm)
{
    int result = MPI_Comm_dup(comm, &communicator);
    Check (result == MPI_SUCCESS);
}

} // end namespace rtt_c4

#endif // C4_MPI

#endif                          // c4_C4_MPI_hh

//---------------------------------------------------------------------------//
//                              end of c4/C4_MPI.hh
//---------------------------------------------------------------------------//
