//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    rng/Random_Inline.hh
 * \author  Paul Henning
 * \brief   Header to bring in appropriate random number generators
 * \note    Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *          All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_rng_Random_Inline_hh
#define rtt_rng_Random_Inline_hh

#include "Rnd_Control_Inline.hh"

namespace rtt_rng {

/*! \brief rn_stream is not used for anything in this library.  It is simply a
 * global variable that some applications use for holding a stream number.
 */
DLL_PUBLIC_rng extern uint64_t rn_stream;
}

#endif // rtt_rng_Random_Inline_hh

//---------------------------------------------------------------------------//
// end of rng/Random_Inline.hh
//---------------------------------------------------------------------------//
