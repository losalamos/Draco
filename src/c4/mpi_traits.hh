//----------------------------------*-C++-*----------------------------------//
// mpi_traits.hh
// Geoffrey Furnish
// Fri Oct  3 09:18:33 1997
//---------------------------------------------------------------------------//
// @> MPI traits specializations.
//---------------------------------------------------------------------------//

#ifndef __c4_mpi_traits_hh__
#define __c4_mpi_traits_hh__

#include "config.hh"

C4_NAMESPACE_BEG

//===========================================================================//
// class mpi_traits - Traits support for the MPI API

// This class and its specializations are intended to support typesafe
// defaulting of the various parameters needed to call MPI API routines.
// Makes life with MPI a lot easier to deal with.
//===========================================================================//

template<class T>
class mpi_traits {
};

template<> class mpi_traits<char> {
  public:
    static MPI_Datatype element_type() { return MPI_CHAR; }
};

template<> class mpi_traits<int> {
  public:
    static MPI_Datatype element_type() { return MPI_INT; }
};

template<> class mpi_traits<float> {
  public:
    static MPI_Datatype element_type() { return MPI_FLOAT; }
};

template<> class mpi_traits<double> {
  public:
    static MPI_Datatype element_type() { return MPI_DOUBLE; }
};

C4_NAMESPACE_END

#endif                          // __c4_mpi_traits_hh__

//---------------------------------------------------------------------------//
//                              end of c4/mpi_traits.hh
//---------------------------------------------------------------------------//
