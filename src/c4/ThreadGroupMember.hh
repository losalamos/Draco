//----------------------------------*-C++-*----------------------------------//
// ThreadGroupMember.hh
// Geoffrey M. Furnish
// Tue Jul 14 11:21:45 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_ThreadGroupMember_hh__
#define __c4_ThreadGroupMember_hh__

#include "config.hh"

#include <pthread.h>

C4_NAMESPACE_BEG

class ThreadControl;

//===========================================================================//
// class ThreadGroupMember - 

// 
//===========================================================================//

class ThreadGroupMember {

  protected:
    pthread_t tid;

    int thread;
    int threads;

    ThreadControl *tcb;

  public:
    ThreadGroupMember();
    virtual ~ThreadGroupMember() {}

    virtual void c4_run_thread() =0;

    template<class X> friend class ThreadGroup;

    void gsync() const;

    template<class T>
    void gsum( T& x ) const;
};

C4_NAMESPACE_END

#include "ThreadGroupMember.t.cc"

#endif                          // __c4_ThreadGroupMember_hh__

//---------------------------------------------------------------------------//
//                              end of c4/ThreadGroupMember.hh
//---------------------------------------------------------------------------//
