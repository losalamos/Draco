//----------------------------------*-C++-*----------------------------------//
// global.cc
// Geoffrey Furnish
// Tue Dec 20 1994
//---------------------------------------------------------------------------//
// @> Global functions provided by the C4 Messaging API.
//---------------------------------------------------------------------------//

// c4 configure
#include <c4/config.h>
// we don't need to compile this file because global.hh includes mpi.t.hh
// if we use KCC
#ifdef C4_NOKCC

#include "global.hh"

// Since global.hh has template function declarations, this file will be
// searched for template function definitions under some compilation
// systems.

#ifdef C4_MPI
#include "mpi.t.hh"
#endif

#endif // C4_NOKCC

//---------------------------------------------------------------------------//
//                              end of global.cc
//---------------------------------------------------------------------------//
