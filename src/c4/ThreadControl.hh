//----------------------------------*-C++-*----------------------------------//
// ThreadControl.hh
// Geoffrey M. Furnish
// Tue Jul 14 11:25:52 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_ThreadControl_hh__
#define __c4_ThreadControl_hh__

#include <pthread.h>

#include "config.hh"

C4_NAMESPACE_BEG

//===========================================================================//
// class ThreadControl - 

// 
//===========================================================================//

class ThreadControl {

// Because we have mutex and condition variable members, we cannot allow
// copying. 

    ThreadControl( const ThreadControl& );
    ThreadControl& operator=( const ThreadControl& );

// Members needed for coordinated actions amongst members of a thread group.

    pthread_mutex_t mutex;
    pthread_cond_t cv;
    int nthreads;
    int cycle;
    int cnt;

    void *scratch;

  public:
    ThreadControl( int nthreads_ );
    ~ThreadControl();

    template<class X> friend class ThreadGroup;
    friend class ThreadGroupMember;
};

C4_NAMESPACE_END

#endif                          // __c4_ThreadControl_hh__

//---------------------------------------------------------------------------//
//                              end of c4/ThreadControl.hh
//---------------------------------------------------------------------------//
