/*!
  \file    rng/Random_Inline.hh
  \author  Paul Henning
  \brief   Header to bring in appropriate random number generators
  \note    Copyright (C) 2008-2011 Los Alamos National Security, LLC.
  \version $Id$
*/

#ifndef rtt_rng_Random_Inline_hh
#define rtt_rng_Random_Inline_hh

//#include "rng/config.h"
#include "Rnd_Control_Inline.hh"

namespace rtt_rng
{

typedef LF_Gen Random;

/*! rn_stream is not used for anything in this library.  It is simply a global
 * variable that some applications use for holding a stream number.
 */
extern int rn_stream;  
}

#endif                          // rtt_rng_Random_Inline_hh

//---------------------------------------------------------------------------//
//                              end of rng/Random_Inline.hh
//---------------------------------------------------------------------------//


