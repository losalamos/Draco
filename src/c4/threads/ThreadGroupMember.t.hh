//----------------------------------*-C++-*----------------------------------//
// ThreadGroupMember.t.hh
// Geoffrey M. Furnish
// Tue Jul 14 11:21:46 1998
//---------------------------------------------------------------------------//
// @> Implements out-of-line member templates of ThreadGroupMember.
//---------------------------------------------------------------------------//

// This file is included directly by ThreadGroupMember.hh, so all the
// includes have already been performed.

#include <stdio.h>

C4_NAMESPACE_BEG

//---------------------------------------------------------------------------//
// Perform a scalar reduction across a family of threads.  The precise
// reduction operation performed is unimportant, and is provided via a
// templating parameter.
//---------------------------------------------------------------------------//

template<class T, class ScalarOp>
void ThreadGroupMember::scalar_reduce( T& x ) const
{
    using namespace std;

// First we need to acquire the mutex.

    int status = pthread_mutex_lock( &tcb->mutex );
    if (status) cerr << "Couldn't lock mutex.\n";

// We've locked the thread control struct, so now we can store the address
// where we need the result stored at the end.

    tcb->v[thread] = &x;

// Now keep a local copy of the cycle index so we can recognize an authentic
// wake-up signal.

    int cycle = tcb->cycle;

    if (tcb->cnt == threads)
    {
    // If we're the first thread, set the global reduction scalar to be our
    // result value.

        tcb->scratch = &x;

    // Also, the first thread needs to initialize the result.

        ScalarOp::init( *static_cast<T *>(tcb->scratch), x );
    }
    else
    {
    // Threads other than the first need to contribute their partial result.

        ScalarOp::apply( *static_cast<T *>(tcb->scratch), x );
    }

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

#ifdef ARRAY_REDUCE_1
template<class T>
void ThreadGroupMember::gsum( T *px, int n ) const
{
    using namespace std;

// First we need to acquire the mutex.

    int status = pthread_mutex_lock( &tcb->mutex );
    if (status) cerr << "Couldn't lock mutex.\n";

// We've locked the thread control struct, so now we can store the address
// where we need the result stored at the end.

    tcb->v[thread] = px;

// Now keep a local copy of the cycle index so we can recognize an authentic
// wake-up signal.

    int cycle = tcb->cycle;

    if (tcb->cnt == threads)
    {
    // If we're the first thread, set the global reduction scalar to be our
    // result value.

        tcb->scratch = px;
    }
    else
    {
    // Threads other than the first need to contribute their partial result.

        for( int i=0; i < n; i++ )
            static_cast<T *>(tcb->scratch)[i] += px[i];
    }

// Decrement the counter and figure out where we are in the pecking order.

    if (--tcb->cnt == 0)
    {
    // We're the last one to arrive, so reset for next cycle.

        tcb->cycle = !tcb->cycle;
        tcb->cnt = tcb->nthreads;

    // Okay, all the data is in, so store it to the result location for each
    // thread.

        for( int i=0; i < threads; i++ )
            for( int j=0; j < n; j++ )
                static_cast<T*>(tcb->v[i])[j] = 
                    static_cast<T *>(tcb->scratch)[j];

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
#endif

#ifdef ARRAY_REDUCE_2
//---------------------------------------------------------------------------//
// The plan with this approach is that we wait for everyone to arrive, then
// we put all the processors to work on summing up independent sections of
// the data and broadcasting it back out to all the various destination
// locations, then they all wait for everyone to be done, then they are all
// released. 
//---------------------------------------------------------------------------//

template<class T>
void ThreadGroupMember::gsum( T *px, int n ) const
{
// Store my address where others can see it.  Cache coherency means threads
// should be able to write to different memory locations without needing to
// lock the mutex.

    tcb->v[thread] = px;

// Now determine the portions of the array that I'm responsible for.  This
// also doesn't need to be inside the mutex.

    int range = n / threads;
    int min = thread * range;
    int max = min + range - 1;
    if (range * threads < n && thread == threads - 1)
        max = n - 1;

// Lock entry mutex.
    int status = pthread_mutex_lock( &tcb->mutex );
    if (status) cerr << "Couldn't lock mutex.\n";

    if (--tcb->cnt == 0)
    {
        Check( tcb->arrival_predicate == 1 );

    // First reset the count.
        tcb->cnt = threads;

    // Also reset the exit predicate.
        tcb->exit_predicate = 1;

    // Now invert the arrival predicate.
        tcb->arrival_predicate = 0;

    // Now wake everyone up.
        status = pthread_cond_broadcast( &tcb->arrival_cv );
        if (status) cerr << "Couldn't broadcast cv.\n";
    }
    else
    {
    // We're not the last one in, so we need to wait till everyone's here.

        while( tcb->arrival_predicate )
        {
            status = pthread_cond_wait( &tcb->arrival_cv, &tcb->mutex );
            if (status) cerr << "Error waiting on cv.\n";
        }
    }

// I've been awoken, and I'm holding the mutex.  Drop it, so someone else can 
// be activated too.

    pthread_mutex_unlock( &tcb->mutex );

// Now sum up the results.

    T *v0 = static_cast<T *>( tcb->v[0] );

    for( int i=1; i < threads; i++ )
    {
        T *vi = static_cast<T *>( tcb->v[i] );
        for( int j=min; j <= max; j++ )
            v0[j] += vi[j];
    }

// Now broadcast the results to all locations.

    for( int i=1; i < threads; i++ )
    {
        T *vi = static_cast<T *>( tcb->v[i] );
        for( int j=min; j <= max; j++ )
            vi[j] = v0[j];
    }

// lock the mutex which indicates I'm finished.  derement the count.  If I'm
// the last one in, wake everybody up, otherwise wait on the signal that says 
// everyone's done.

    pthread_mutex_lock( &tcb->mutex );

    if (--tcb->cnt == 0)
    {
        Check( tcb->exit_predicate == 1 );

    // First reset the count.
        tcb->cnt = threads;

    // Also reset the arrival_predicate.
        tcb->arrival_predicate = 1;

    // Now flip the exit predicate.
        tcb->exit_predicate = 0;

    // Now wake everyone up.
        status = pthread_cond_broadcast( &tcb->exit_cv );
        if (status) cerr << "Couldn't broadcast cv.\n";
    }
    else
    {
    // Have to wait on everyone to finish reducing.

        while( tcb->exit_predicate )
        {
            status = pthread_cond_wait( &tcb->exit_cv, &tcb->mutex );
            if (status) cerr << "Error waiting on cv.\n";
        }
    }

// I've been awoken, and I'm holding the mutex.  All the work is done, so
// just drop the mutex and skeedaddle.

    pthread_mutex_unlock( &tcb->mutex );
}
#endif

#ifdef ARRAY_REDUCE_3
//---------------------------------------------------------------------------//
// This one is based on ARRAY_REDUCE_2, but seeks to overcome the cache
// coherency bug which was discovered in that one.  We do this by having each 
// thread store its results to its subregion of its own buffer (instead of
// every thread summing to its subregion of thread 0's buffer).  Then we have 
// to introduce a synchronization point to guarantee that all the threads see 
// the new values, then every thread loads its buffer by fetching data from
// the subregions of the other thread buffers that they stored into.  Sheesh.
//---------------------------------------------------------------------------//

template<class T>
void ThreadGroupMember::gsum( T *px, int n ) const
{
    Require( tcb->arrival_predicate == 1 );
    Require( tcb->middle_predicate == 1 );

// First let's calculate the region of the buffer that we're responsible
// for.  We don't need the mutex to calculate this.

    int range = n / threads;
    int min = thread * range;
    int max = min + range - 1;
    if (range * threads < n && thread == threads - 1)
        max = n - 1;

// Lock the entry mutex.

    int status = pthread_mutex_lock( &tcb->mutex );
    if (status) cerr << "Couldn't lock mutex.\n";

// Store our data into tcb where everyone can see it.  Note that each thread
// is storing data in different memory locations, so you ought to be able to
// do this outside the mutex.  But, because of the SGI cache coherency bug
// discovered in ARRAY_REDUCE_2 above, we have to guard against storing data
// into tcb->v[i] unless protected by a mutex.  This RRRRRRRRREALLY sucks.

    tcb->v[thread] = px;
    tcb->min[thread] = min;
    tcb->max[thread] = max;

    int arrival_predicate = tcb->arrival_predicate;

    if (--tcb->cnt == 0)
    {
        Check( tcb->arrival_predicate == 1 );

    // Configure the exit predicate.
        tcb->exit_predicate = 1;

    // Now reset the count.
        tcb->cnt = threads;

    // Now invert the arrival predicate.
        tcb->arrival_predicate = 0;

    // Now wake everyone up.
        status = pthread_cond_broadcast( &tcb->arrival_cv );
        if (status) cerr << "Couldn't broadcast cv.\n";
    }
    else
    {
    // We're not the last one in, so we need to wait till everyone's here.

        while( arrival_predicate == tcb->arrival_predicate )
        {
            status = pthread_cond_wait( &tcb->arrival_cv, &tcb->mutex );
            if (status) cerr << "Error waiting on cv.\n";
        }
    }

// I've been awoken, and I'm holding the mutex.  Drop it, so someone else can 
// be activated too.

    pthread_mutex_unlock( &tcb->mutex );

// Now sum up the results.

    for( int i=0; i < threads; i++ )
    {
        if (i == thread) continue;
        T *vi = static_cast<T *>( tcb->v[i] );
        for( int j=min; j <= max; j++ )
            px[j] += vi[j];
    }

// Now we need to have everyone wait on a condition variable in order to
// guarantee that the results calculated by each thread are visible by all.

    pthread_mutex_lock( &tcb->mutex );

    int middle_predicate = tcb->middle_predicate;

    if ( --tcb->cnt == 0 )
    {
        Check( tcb->middle_predicate == 1 );

    // Reset the arrival predicate so it will be valid on entry next time.
        tcb->arrival_predicate = 1;

    // Reset the count so the next wait works.
        tcb->cnt = threads;

    // Now flip the predicate for this wait.
        tcb->middle_predicate = 0;

    // Now wake everyone up.
        status = pthread_cond_broadcast( &tcb->middle_cv );
        if (status) cerr << "Couldn't broadcast cv.\n";
    }
    else
    {
    // We're not the last one in, so we need to wait till everyone's here.

        while( middle_predicate == tcb->middle_predicate )
        {
            status = pthread_cond_wait( &tcb->middle_cv, &tcb->mutex );
            if (status) cerr << "Error waiting on cv.\n";
        }
    }

// We're through the middle sync point, and we have the mutex.  Unlock it so
// others can proceed.

    pthread_mutex_unlock( &tcb->mutex );

// Now go fetch the results from each thread.

    for( int i=0; i < threads; i++ )
    {
        if (i == thread) continue;
        T *vi = static_cast<T *>( tcb->v[i] );
        for( int j=tcb->min[i]; j <= tcb->max[i]; j++ )
            px[j] = vi[j];
    }

// Lock the exit mutex and register termination condition.

    pthread_mutex_lock( &tcb->mutex );

    int exit_predicate = tcb->exit_predicate;

    if (--tcb->cnt == 0)
    {
        Check( tcb->exit_predicate == 1 );

    // Make sure the entry predicate is set back to starting value.

        tcb->middle_predicate = 1;

    // First reset the count.

        tcb->cnt = threads;

    // Now flip the exit predicate.

        tcb->exit_predicate = 0;

    // Now wake everyone up.

        status = pthread_cond_broadcast( &tcb->exit_cv );
        if (status) cerr << "Couldn't broadcast cv.\n";
    }
    else
    {
    // Have to wait on everyone to finish fetching the results.

        while( exit_predicate == tcb->exit_predicate )
        {
            status = pthread_cond_wait( &tcb->exit_cv, &tcb->mutex );
            if (status) cerr << "Error waiting on cv.\n";
        }
    }

// I've been awoken, and I'm holding the mutex.  All the work is done, so
// just drop the mutex and skeedaddle.

    pthread_mutex_unlock( &tcb->mutex );
}
#endif

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of ThreadGroupMember.t.hh
//---------------------------------------------------------------------------//
