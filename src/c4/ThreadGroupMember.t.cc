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

template<class T>
void ThreadGroupMember::gsum( T& x ) const
{
    T local_sum = T();          // #'s init to zero this way.

// First we need to acquire the mutex.
    int status = pthread_mutex_lock( &tcb->mutex );
    if (status) cerr << "Couldn't lock mutex.\n";

// We've locked the thread control struct, so now we can store the address
// where we need the result stored at the end.

    tcb->v[thread] = &x;

// Now keep a local copy of the cycle index so we can recognize an authentic
// wake-up signal.
    int cycle = tcb->cycle;

// If we're the first thread, set the scratch pointer to point to our
// function scope local_sum variable.  :-).

    if (tcb->cnt == threads)
        tcb->scratch = &local_sum;

// Whichever thread we are, we can now go ahead and add our partial sum into
// the total.

    *static_cast<T *>(tcb->scratch) += x;

// Decrement the counter and figure out where we are in the pecking order.
    if (--tcb->cnt == 0)
    {
    // We're the last one to arrive, so reset for next cycle.

        tcb->cycle = !tcb->cycle;
        tcb->cnt = tcb->nthreads;

    // Okay, all the data is in, so store it to the result location for each
    // thread.

        for( int i=0; i < threads; i++ )
            *static_cast<T*>(tcb->v[i]) = 
                *static_cast<T *>(tcb->scratch);

    // Now we can wake everybody up.
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

    // Note that the waking thread has already propagated the result value to 
    // our parameter x, so there is nothing left for us to do.

    // In particular, note that even though we now hold the mutex on the
    // thread control block tcb, it would not be valid to pull some result
    // out of tcb at this point since some other member of our thread group,
    // woken beforehand, may have already reached its next reduction point,
    // and may have already invalidated the data in tcb.  So, given this
    // wakeup strategy, it is really essential that /all/ the work, already
    // be done at this point.
    }

// Everybody wakes up with the mutex locked, so unlock it now.
    status = pthread_mutex_unlock( &tcb->mutex );
    if (status) cerr << "error unlocking mutex.\n";
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of ThreadGroupMember.cc
//---------------------------------------------------------------------------//
