//----------------------------------*-C++-*----------------------------------//
// global_scalar.hh
// Geoffrey Furnish
// 12 November 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_global_scalar_hh__
#define __c4_global_scalar_hh__

#include "ds++/Assert.hh"

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Blocking send and receive functions.
//---------------------------------------------------------------------------//

template<class T>
int Send( const T *buf, int nels, int dest,
	  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    return C4_SUCCESS;
}

template<class T>
int Recv( T *buf, int nels, int source,
	  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    return 0;
}

//---------------------------------------------------------------------------//
// Asynchronous send and recive functions, returning a request handle.
//---------------------------------------------------------------------------//

template<class T>
C4_Req SendAsync( const T *buf, int nels, int dest,
		  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    C4_Req r;
    return r;
}

template<class T>
C4_Req RecvAsync( T *buf, int nels, int source,
		  int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    C4_Req r;
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
}

template<class T>
void RecvAsync( C4_Req& r, T *buf, int nels, int source,
		int tag /*=c4_traits<T*>::Tag*/, int group /*=0*/ )
{
    Require( !r.inuse() );
}

template<class T> void gsum( T& x ) {}
template<class T> void gmin( T& x ) {}
template<class T> void gmax( T& x ) {}

template<class T> void gsum( T *px, int n ) {}
template<class T> void gmin( T *px, int n ) {}
template<class T> void gmax( T *px, int n ) {}

C4_NAMESPACE_END

#endif                          // __c4_global_scalar_hh__

//---------------------------------------------------------------------------//
//                              end of c4/global_scalar.hh
//---------------------------------------------------------------------------//
