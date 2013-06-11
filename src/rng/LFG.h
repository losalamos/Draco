/*----------------------------------*-C-*----------------------------------*/
/*!
 * \file   rng/LFG.h
 * \author Paul Henning
 * \date   June 28, 2006
 * \brief  Lagged Fibonacci Generator Random Number Generator
 * \note   Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
/*---------------------------------------------------------------------------*/
/* $Id$                       */
/*---------------------------------------------------------------------------*/

#ifndef rtt_rng_LFG_H
#define rtt_rng_LFG_H

#include "rng/config.h"
#include "ds++/config.h"

#ifdef __cplusplus
extern "C" { 
#endif
    
#ifndef LFG_PARAM_SET
#  ifdef RNG_NR
#    define LFG_PARAM_SET 2
#  else 
#    define LFG_PARAM_SET 1
#  endif
#endif 
    
/* 
   This is kludgy.  LFG_DATA_SIZE is the number of unsigned ints required to
   store the state of a random number generator. This needs to be the same as
   the value "OFFSET_END" in LFG.c.  
*/

#if LFG_PARAM_SET == 1
#  define LFG_DATA_SIZE 54
#elif LFG_PARAM_SET == 2
#  define LFG_DATA_SIZE 96
#endif

/* 
   Create an new random number generator.  This is generator "gennum", which
   is in the range [0, total_gen).  The seed value should be constant for all
   calls to lfg_create_rng.  The random number state is written into the
   memory addressed by [*begin, *end), where (end-begin > LFG_DATA_SIZE), and
   you are responsible for managing that memory (both allocation and free).

   NOTE: [begin, end) MUST be contiguous addresses!
*/
extern void lfg_create_rng(const unsigned gennum, 
                           unsigned seed,
                           unsigned* begin,
                           unsigned* end);


/* 
   Same as above, but doesn't complete the initialization.  You MUST call
   part two!

   NOTE: [begin, end) MUST be contiguous addresses!
*/
extern void lfg_create_rng_part1(const unsigned gennum, 
                                 unsigned seed,
                                 unsigned* begin,
                                 unsigned* end);

extern void lfg_create_rng_part2(unsigned* begin
                                 /*, unsigned* end */);

/*
  Create a new random number generator from an already existing one.  This
  produces an independent stream.  Data is written in to [*begin, *end), see
  lfg_create_rng() for more details.
*/
extern void lfg_spawn_rng(unsigned* genptr, 
			  unsigned* begin, 
			  unsigned* end);

/*
  Get the next double from the random number generator pointed to by genptr
  (this was probably "begin" in one of the creation functions).
*/
extern double lfg_gen_dbl(unsigned* genptr);

/*
  Get the next integer from the random number generator pointed to by genptr
  (this was probably "begin" in one of the creation functions).  
*/
extern int lfg_gen_int(unsigned* genptr);

/*
  Dump a little diagnostics printout.
*/
extern void lfg_print(unsigned* genptr);

/*
  Return the number of unsigned ints needed to hold the state.  This should 
  be the same as the LFG_DATA_SIZE #define. 
*/
extern unsigned lfg_size();

/*
  Return the gennum associated with this stream.
*/
extern unsigned lfg_gennum(unsigned* genptr);

/*
 Return a "unique" number associated with this stream and state
*/
extern unsigned lfg_unique_num(unsigned* genptr);

extern DLL_PUBLIC 
void errprint(char const *level, char const *routine, char const *error);

#ifdef __cplusplus
}
#endif

#endif
