//----------------------------------*-C++-*----------------------------------//
// tags.hh
// Geoffrey Furnish
// Fri Oct  3 10:06:44 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_tags_hh__
#define __c4_tags_hh__

#include "c4/config.hh"

const int C4_SUCCESS = 0;
#ifdef __PARAGON__
const int C4_Any_Tag = -1;
const int C4_Any_Source = -1;
#endif
#ifdef __C4_SCALAR__
const int C4_Any_Tag = -1;
const int C4_Any_Source = -1;
#endif
#ifdef C4_SHMEM
const int C4_Any_Tag = -1;
const int C4_Any_Source = -1;
#endif
#ifdef __MPI__
const int C4_Any_Tag = MPI_ANY_TAG;
const int C4_Any_Source = MPI_ANY_SOURCE;
#endif

const int C4_int_Tag = 432;
const int C4_float_Tag = 433;
const int C4_double_Tag = 434;
const int C4_int_ptr_Tag = 443;
const int C4_float_ptr_Tag = 444;
const int C4_double_ptr_Tag = 445;

enum { C4_Pass_Left, C4_Pass_Right };

#endif                          // __c4_tags_hh__

//---------------------------------------------------------------------------//
//                              end of c4/tags.hh
//---------------------------------------------------------------------------//
