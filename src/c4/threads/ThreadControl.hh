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

#include "../config.hh"

#include <vector>

// #define ARRAY_REDUCE_1
#define ARRAY_REDUCE_2
// #define ARRAY_REDUCE_3

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

    std::vector<void *> v;

#ifdef ARRAY_REDUCE_2
    int arrival_predicate;
    pthread_cond_t arrival_cv;
    int exit_predicate;
    pthread_cond_t exit_cv;
#endif

#ifdef ARRAY_REDUCE_3
    int arrival_predicate;
    pthread_cond_t arrival_cv;

    int middle_predicate;
    pthread_cond_t middle_cv;

    int exit_predicate;
    pthread_cond_t exit_cv;

    std::vector<int> min, max;
#endif

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
