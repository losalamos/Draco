//----------------------------------*-C++-*----------------------------------//
// mpi_traits.hh
// Geoffrey Furnish
// Fri Oct  3 09:18:33 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_mpi_traits_hh__
#define __c4_mpi_traits_hh__

//===========================================================================//
// class mpi_traits - 

// 
//===========================================================================//

template<class T>
class mpi_traits {
};

template<> class mpi_traits<int *> {
  public:
    static const MPI_Datatype element_type = MPI_INT;
};

template<> class mpi_traits<float *> {
  public:
    static const MPI_Datatype element_type = MPI_FLOAT;
};

template<> class mpi_traits<double *> {
  public:
    static const MPI_Datatype element_type = MPI_DOUBLE;
};

#endif                          // __c4_mpi_traits_hh__

//---------------------------------------------------------------------------//
//                              end of c4/mpi_traits.hh
//---------------------------------------------------------------------------//
