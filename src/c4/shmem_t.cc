//----------------------------------*-C++-*----------------------------------//
// shmem_t.cc
// Geoffrey M. Furnish
// Wed Feb 18 15:55:22 1998
//---------------------------------------------------------------------------//
// @> Template functions for the C4 MPI interface.
//---------------------------------------------------------------------------//

#include "c4_traits.hh"

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
    return SHM_Send( (void *) buf, nels * sizeof(T), dest, tag, group );
}

template<class T>
int Recv( T *buf, int nels, int source,
	  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    int bytes = SHM_Recv( (void *) buf, nels*sizeof(T), source, tag, group );
    return bytes/sizeof(T);
}

//---------------------------------------------------------------------------//
// Asynchronous send and recive functions, returning a request handle.
//---------------------------------------------------------------------------//

template<class T>
C4_Req SendAsync( const T *buf, int nels, int dest,
		  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    return SHM_SendAsync( (void *) buf, nels*sizeof(T), dest, tag, group );
}

template<class T>
C4_Req RecvAsync( T *buf, int nels, int source,
		  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    return SHM_RecvAsync( (void *) buf, nels*sizeof(T), source, tag, group );
}

//---------------------------------------------------------------------------//
// Asynchronous send and recive functions, using a provided request handle.
//---------------------------------------------------------------------------//

template<class T>
void SendAsync( C4_Req& r, const T *buf, int nels, int dest,
		int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    Require( !r.inuse() );
//     r.set();
//     MPI_Isend( (void *) buf, nels, mpi_traits<T>::element_type>,
// 	       dest, tag, MPI_COMM_WORLD, &r.r );
}

template<class T>
void RecvAsync( C4_Req& r, T *buf, int nels, int source,
		int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    Require( !r.inuse() );
//     r.set();
//     MPI_Irecv( (void *) buf, nels, mpi_traits<T>::element_type>,
// 	       source, tag, MPI_COMM_WORLD, &r.r );
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of shmem_t.cc
//---------------------------------------------------------------------------//
