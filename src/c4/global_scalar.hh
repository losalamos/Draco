//----------------------------------*-C++-*----------------------------------//
// global_scalar.hh
// Geoffrey Furnish
// 12 November 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_global_scalar_hh__
#define __c4_global_scalar_hh__

C4_NAMESPACE_BEG

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

C4_NAMESPACE_END

#endif                          // __c4_global_scalar_hh__

//---------------------------------------------------------------------------//
//                              end of c4/global_scalar.hh
//---------------------------------------------------------------------------//
