//----------------------------------*-C++-*----------------------------------//
// ThreadControl.cc
// Geoffrey M. Furnish
// Tue Jul 14 11:25:53 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/ThreadControl.hh"

#include <iostream>

C4_NAMESPACE_BEG

ThreadControl::ThreadControl( int nthreads_ )
    : nthreads( nthreads_ ),
      cycle( 0 ),
      cnt( nthreads ),
      v( nthreads )
{
    if (pthread_mutex_init( &mutex, NULL ))
        throw "Can't initialize mutex.";

    if (pthread_cond_init( &cv, NULL ))
        throw "Can't initialize condition variable.";
}

ThreadControl::~ThreadControl()
{
    if (pthread_mutex_destroy( &mutex ))
        throw "Can't destroy mutex.";

    if (pthread_cond_destroy( &cv ))
        throw "Can't destroy condition variable.";
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of ThreadControl.cc
//---------------------------------------------------------------------------//
