//----------------------------------*-C++-*----------------------------------//
// global_mpi.hh
// Geoffrey Furnish
// 12 November 1997
//---------------------------------------------------------------------------//
// @> C4 over MPI header.
//---------------------------------------------------------------------------//

#ifndef __c4_global_mpi_hh__
#define __c4_global_mpi_hh__

namespace C4 {

    inline double Wtime() { return MPI_Wtime(); }
    inline double Wtick() { return MPI_Wtick(); }

}

#include "mpi.t.hh"

#endif                          // __c4_global_mpi_hh__

//---------------------------------------------------------------------------//
//                              end of c4/global_mpi.hh
//---------------------------------------------------------------------------//
