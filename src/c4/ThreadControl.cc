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
#ifdef ARRAY_REDUCE_3
    , min( nthreads ), max (nthreads )
#endif
{
    if (pthread_mutex_init( &mutex, NULL ))
        throw "Can't initialize mutex.";

    if (pthread_cond_init( &cv, NULL ))
        throw "Can't initialize condition variable.";

#ifdef ARRAY_REDUCE_2
    arrival_predicate = 1;

    if (pthread_cond_init( &arrival_cv, NULL ))
        throw "Can't initialize condition variable.";

    exit_predicate = 1;

    if (pthread_cond_init( &exit_cv, NULL ))
        throw "Can't initialize condition variable.";
#endif

#ifdef ARRAY_REDUCE_3
    arrival_predicate = 1;

    if (pthread_cond_init( &arrival_cv, NULL ))
        throw "Can't initialize condition variable.";

    middle_predicate = 1;

    if (pthread_cond_init( &middle_cv, NULL ))
        throw "Can't initialize condition variable.";

    exit_predicate = 1;

    if (pthread_cond_init( &exit_cv, NULL ))
        throw "Can't initialize condition variable.";
#endif
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
