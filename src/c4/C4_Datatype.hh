//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/C4_Datatype.hh
 * \author Kent G. Budge
 * \brief  C4_Datatype class definition.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef c4_C4_Datatype_hh
#define c4_C4_Datatype_hh

// C4 package configure
#include "c4/config.h"

#ifdef C4_MPI
#include "c4_mpi.h"
#endif

namespace rtt_c4 {

#ifdef C4_MPI

typedef MPI_Datatype C4_Datatype;

#else

// If serial, make this a brain-dead type. It won't actually be used.

typedef void *C4_Datatype;

#endif

} // end namespace rtt_c4

#endif // c4_C4_Datatype_hh

//----------------------------------------------------------------------------//
// end of c4/C4_Datatype.hh
//----------------------------------------------------------------------------//
