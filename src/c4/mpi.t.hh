//----------------------------------*-C++-*----------------------------------//
// mpi_t.cc
// Geoffrey Furnish
// Fri Oct  3 09:58:34 1997
//---------------------------------------------------------------------------//
// @> Template functions for the C4 MPI interface.
//---------------------------------------------------------------------------//

#include "mpi_traits.hh"
#include "c4_traits.hh"
#include "ds++/Assert.hh"

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// MPI send/receive calls (basic set)
//
// Synchronous:
//	MPI_Send(void* buf, int count, MPI_Datatype datatype,
//		 int dest, int tag, MPI_Comm comm);
//	MPI_Recv(void* buf, int count, MPI_Datatype datatype,
//		 int source, int tag, MPI_Comm comm, MPI_Status *status);
//
// Asynchronous:
//	MPI_Isend(void* buf, int count, MPI_Datatype datatype,
//		  int dest, int tag, MPI_Comm comm, MPI_Request *request);
//	MPI_Irecv(void* buf, int count, MPI_Datatype datatype,
//		  int source, int tag, MPI_Comm comm, MPI_Request *request);
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Basic send and receive functions.
//---------------------------------------------------------------------------//

template<class T>
int Send( const T *buf, int nels, int dest,
	  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    MPI_Send( const_cast<T *>(buf), nels, mpi_traits<T>::element_type(),
	      dest, tag, MPI_COMM_WORLD );
    return C4_SUCCESS;
}

template<class T>
int Recv( T *buf, int nels, int source,
	  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    int cnt;
    MPI_Status status;
    MPI_Recv( buf, nels, mpi_traits<T>::element_type(),
	      source, tag, MPI_COMM_WORLD, &status );
    MPI_Get_count( &status, mpi_traits<T>::element_type(), &cnt );
    return cnt;
}

//---------------------------------------------------------------------------//
// Asynchronous send and recive functions, returning a request handle.
//---------------------------------------------------------------------------//

template<class T>
C4_Req SendAsync( const T *buf, int nels, int dest,
		  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    C4_Req r;
    MPI_Isend( (void *) buf, nels, mpi_traits<T>::element_type(),
	       dest, tag, MPI_COMM_WORLD, &r.r() );
    r.set();
    return r;
}

template<class T>
C4_Req RecvAsync( T *buf, int nels, int source,
		  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    C4_Req r;
    MPI_Irecv( (void *) buf, nels, mpi_traits<T>::element_type(),
	       source, tag, MPI_COMM_WORLD, &r.r() );
    r.set();
    return r;
}

//---------------------------------------------------------------------------//
// Asynchronous send and recive functions, using a provided request handle.
//---------------------------------------------------------------------------//

template<class T>
void SendAsync( C4_Req& r, const T *buf, int nels, int dest,
		int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    Require( !r.inuse() );
    r.set();
    MPI_Isend( (void *) buf, nels, mpi_traits<T>::element_type(),
	       dest, tag, MPI_COMM_WORLD, &r.r() );
}

template<class T>
void RecvAsync( C4_Req& r, T *buf, int nels, int source,
		int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    Require( !r.inuse() );
    r.set();
    MPI_Irecv( (void *) buf, nels, mpi_traits<T>::element_type(),
	       source, tag, MPI_COMM_WORLD, &r.r() );
}

template<class T>
void gsum( T& x )
{
    T y = x;
    MPI_Allreduce( &y, &x, 1, mpi_traits<T>::element_type(),
		   MPI_SUM, MPI_COMM_WORLD );
}

template<class T>
void gprod( T& x )
{
    T y = x;
    MPI_Allreduce( &y, &x, 1, mpi_traits<T>::element_type(),
		   MPI_PROD, MPI_COMM_WORLD );
}

template<class T>
void gmin( T& x )
{
    T y = x;
    MPI_Allreduce( &y, &x, 1, mpi_traits<T>::element_type(),
		   MPI_MIN, MPI_COMM_WORLD );
}

template<class T>
void gmax( T& x )
{
    T y = x;
    MPI_Allreduce( &y, &x, 1, mpi_traits<T>::element_type(),
		   MPI_MAX, MPI_COMM_WORLD );
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of mpi_t.cc
//---------------------------------------------------------------------------//
