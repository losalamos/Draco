//----------------------------------*-C++-*----------------------------------//
// global.hh
// Geoffrey Furnish
// Tue Dec 20 1994
//---------------------------------------------------------------------------//
// @> Global functions provided by the C4 Messaging API.
//---------------------------------------------------------------------------//

#ifndef __c4_global_hh__
#define __c4_global_hh__

// c4 package configure
#include <c4/config.h>

#include "config.hh"
#include "tags.hh"
#include "c4_traits.hh"

namespace C4
{

// Forward Reference

class C4_Req;

void Init( int& argc, char **& argv );
void Finalize();

// Informational 

int node();
int nodes();
int group();

// Global sync

void gsync();

// Send/receive

template<class T>
int Send( const T *buf, int nels, int dest,
	  int tag =c4_traits<T*>::Tag, int group =0 );
template<class T>
int Recv( T *buf, int nels, int source,
	  int tag =c4_traits<T*>::Tag, int group =0 );

// Async send/receive

template<class T>
C4_Req SendAsync( const T *buf, int nels, int dest,
		  int tag =c4_traits<T*>::Tag, int group =0 );
template<class T>
C4_Req RecvAsync( T *buf, int nels, int source,
		  int tag =c4_traits<T*>::Tag, int group =0 );

// Convenience forms which avoid object creation/copy costs from handling the 
// C4_Req's.

template<class T>
void SendAsync( C4_Req& r, const T *buf, int nels, int dest,
		int tag =c4_traits<T*>::Tag, int group =0 );
template<class T>
void RecvAsync( C4_Req& r, T *buf, int nels, int source,
		int tag =c4_traits<T*>::Tag, int group =0 );

// Convenience forms for the elemental send/receive ops.

template<class T>
inline int Send( const T& data, int dest,
		 int tag =c4_traits<T>::Tag, int group =0 )
{
    return Send( &data, 1, dest, tag, group );
}

template<class T>
inline int Recv( T& data, int source,
		 int tag =c4_traits<T>::Tag, int group =0 )
{
    return Recv( &data, 1, source, tag, group );
}

// Global reductions.

// Scalar.

template<class T> void gsum( T& x );
template<class T> void gprod( T& x );
template<class T> void gmin( T& x );
template<class T> void gmax( T& x );

// Array.

template<class T> void gsum( T *px, int n );
template<class T> void gprod( T *px, int n );
template<class T> void gmin( T *px, int n );
template<class T> void gmax( T *px, int n );

// Timing

double Wtime();
double Wtick();

} // end namespace C4

//---------------------------------------------------------------------------//
// Now include anything which might help a particular hardware abstraction
// layer.  For example, C4_gsync might be inlined...

#ifdef C4_SHMEM
#include "global_shmem.hh"
#endif
#ifdef C4_MPI
#include "global_mpi.hh"
#endif
#ifdef C4_SCALAR
#include "global_scalar.hh"
#endif

// For convenience let's include the C4_Req header

#include "C4_Req.hh"

#endif                          // __c4_global_hh__

//---------------------------------------------------------------------------//
//                              end of c4/global.hh
//---------------------------------------------------------------------------//
