//----------------------------------*-C++-*----------------------------------//
// global_mpi.cc
// Maurice LeBrun
// Wed Feb  1 16:50:59 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for MPI
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "c4/mpi_traits.hh"

#include "ds++/DynArray.hh"
#include "ds++/Assert.hh"

using dsxx::DynArray;

//---------------------------------------------------------------------------//
// Miscellaneous

C4_NAMESPACE_BEG

void Init( int& argc, char **& argv )
{
    MPI_Init( &argc, &argv );
}

void Finalize()
{
    MPI_Finalize();
}

int node()
{
    int node;
    MPI_Comm_rank( MPI_COMM_WORLD, &node );
    return node;
}

int nodes()
{
    int nodes;
    MPI_Comm_size( MPI_COMM_WORLD, &nodes );
    return nodes;
}

int group()
{
    int group = 0;
    return group;
}

void gsync()
{
    MPI_Barrier( MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Global reduction operations.
//
// The call most like that in NX is:
//
//	MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
//		      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
// 
// which returns the result to all processes in the group.  "op" determines
// the type of reduction performed.
//
// Available reduction operators, description, and allowed types: 
//
// MPI_MAX, MPI_MIN		min, max			I, F
// MPI_SUM, MPI_PROD		sum, product			I, F
// MPI_BAND, MPI_BOR, MPI_BXOR	bitwise and, or, xor		I, B
// MPI_LAND, MPI_LOR, MPI_LXOR	logical and, or, xor		I
// MPI_MAXLOC, MPI_MINLOC	min, max value and location
//
// where types are:
//
// I:	MPI_INT, MPI_LONG, MPI_SHORT, MPI_UNSIGNED_SHORT, 
//	MPI_UNSIGNED, MPI_UNSIGNED_LONG
// F:	MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE
// B:	MPI_BYTE
//---------------------------------------------------------------------------//

// MPI array reductions traits class.

template<class T> class mpi_ar_traits
{
  public:
    static DynArray<T> ar_buf;
};

template<class T>
void gsum( T *px, int n )
{
    Assert( n >= 0 );
    T dummy = T();

    mpi_ar_traits<T>::ar_buf[n-1] = dummy; // auto expand the buffer.
    for( int i=0; i < n; i++ )
	mpi_ar_traits<T>::ar_buf[i] = px[i];

    MPI_Allreduce( &mpi_ar_traits<T>::ar_buf[0], px, n,
		   mpi_traits<T>::element_type(), MPI_SUM, MPI_COMM_WORLD );
}

template<class T>
void gmin( T *px, int n )
{
    Assert( n >= 0 );
    T dummy = T();

    mpi_ar_traits<T>::ar_buf[n-1] = dummy; // auto expand the buffer.
    for( int i=0; i < n; i++ )
	mpi_ar_traits<T>::ar_buf[i] = px[i];

    MPI_Allreduce( &mpi_ar_traits<T>::ar_buf[0], px, n,
		   mpi_traits<T>::element_type(), MPI_MIN, MPI_COMM_WORLD );
}

template<class T>
void gmax( T *px, int n )
{
    Assert( n >= 0 );
    T dummy = T();

    mpi_ar_traits<T>::ar_buf[n-1] = dummy; // auto expand the buffer.
    for( int i=0; i < n; i++ )
	mpi_ar_traits<T>::ar_buf[i] = px[i];

    MPI_Allreduce( &mpi_ar_traits<T>::ar_buf[0], px, n,
		   mpi_traits<T>::element_type(), MPI_MAX, MPI_COMM_WORLD );
}

template<class T>
DynArray<T> mpi_ar_traits<T>::ar_buf(10);

template void gsum( int *px, int n );
template void gsum( float *px, int n );
template void gsum( double *px, int n );

template void gmin( int *px, int n );
template void gmin( float *px, int n );
template void gmin( double *px, int n );

template void gmax( int *px, int n );
template void gmax( float *px, int n );
template void gmax( double *px, int n );

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of global_mpi.cc
//---------------------------------------------------------------------------//
