//----------------------------------*-C++-*----------------------------------//
// global.cc
// Geoffrey Furnish
// Tue Dec 20 1994
//---------------------------------------------------------------------------//
// @> Global functions provided by the C4 Messaging API.
//---------------------------------------------------------------------------//

#include "c4/global.hh"

// Since global.hh has template function declarations, this file will be
// searched for template function definitions under some compilation
// systems.

#ifdef C4_MPI
#include "c4/mpi_t.cc"
#endif

//---------------------------------------------------------------------------//
//                              end of global.cc
//---------------------------------------------------------------------------//
