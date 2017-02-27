//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/global.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 10:56:16 2002
 * \brief  C4 function declarations and class * definitions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This file allows the client to include the message passing services
 * provided by C4.  The function declarations and class definitions are
 * contained in the rtt_c4 namespace.  For backwards compatibility, the
 * old-style C4 functions and classes are declared in the C4 namespace.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_c4_global_hh
#define rtt_c4_global_hh

// C4 Message Passing Functions
#include "C4_Functions.hh"

//---------------------------------------------------------------------------//
// Include the appropriate header for an underlying message passing
// implementation.  This allows the definition of inline functions declared
// in C4_Functions.hh.

#ifdef C4_SCALAR
#include "C4_Serial.hh"
#endif

#ifdef C4_MPI
#include "C4_MPI.hh"
#endif

#endif // rtt_c4_global_hh

//---------------------------------------------------------------------------//
// end of c4/global.hh
//---------------------------------------------------------------------------//
