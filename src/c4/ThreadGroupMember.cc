//----------------------------------*-C++-*----------------------------------//
// ThreadGroupMember.cc
// Geoffrey M. Furnish
// Tue Jul 14 11:21:46 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/ThreadGroupMember.hh"
#include "c4/ThreadControl.hh"

#include <iostream>
using namespace std;

C4_NAMESPACE_BEG

ThreadGroupMember::ThreadGroupMember() 
{
    tcb = 0;
}

void ThreadGroupMember::gsync() const
{
// First we need to acquire the mutex.
    int status = pthread_mutex_lock( &tcb->mutex );
    if (status) cerr << "Couldn't lock mutex.\n";

    int cycle = tcb->cycle;

    if (--tcb->cnt == 0)
    {
    // We're the last one to arrive, so reset for next cycle, and then wake
    // everybody up.

        tcb->cycle = !tcb->cycle;
        tcb->cnt = tcb->nthreads;

        status = pthread_cond_broadcast( &tcb->cv );
        if (status) cerr << "Couldn't broadcast cv.\n";
    }
    else
    {
    // We're not the last one to arrive, so wait for the predicate to
    // change.

        while( cycle == tcb->cycle )
        {
            status = pthread_cond_wait( &tcb->cv, &tcb->mutex );
            if (status) cerr << "Error waiting on cv.\n";
        }
    }

// Everybody wakes up with the mutex locked, so unlock it now.
    status = pthread_mutex_unlock( &tcb->mutex );
    if (status) cerr << "error unlocking mutex.\n";
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of ThreadGroupMember.cc
//---------------------------------------------------------------------------//
