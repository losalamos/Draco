//----------------------------------*-C++-*----------------------------------//
// mpi_t.cc
// Geoffrey Furnish
// Fri Oct  3 09:58:34 1997
//---------------------------------------------------------------------------//
// @> Template functions for the C4 MPI interface.
//---------------------------------------------------------------------------//

#include "mpi_traits.hh"
#include "c4_traits.hh"

C4_NAMESPACE_BEG

template<class T>
int Send( const T *buf, int nels, int dest, int group /*=0*/ )
{
    MPI_Send( (void *) buf, nels, mpi_traits<T*>::element_type,
	      dest, c4_traits<T*>::Tag, MPI_COMM_WORLD );
    return C4_SUCCESS;
}

template<class T>
int Recv( T *buf, int nels, int source, int group /*=0*/ )
{
    int cnt;
    MPI_Status status;
    MPI_Recv( buf, nels, mpi_traits<T*>::element_type,
	      source, c4_traits<T*>::Tag, MPI_COMM_WORLD, &status );
    MPI_Get_count( &status, mpi_traits<T*>::element_type, &cnt );
    return cnt;
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of mpi_t.cc
//---------------------------------------------------------------------------//
