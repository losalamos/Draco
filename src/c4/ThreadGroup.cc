//----------------------------------*-C++-*----------------------------------//
// ThreadGroup.cc
// Geoffrey M. Furnish
// Mon Jul 13 12:56:20 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/ThreadGroup.hh"
#include "c4/ThreadGroupMember.hh"
#include "c4/ThreadControl.hh"

C4_NAMESPACE_BEG

extern "C" void *thread_start_routine( void *arg )
{
    ThreadGroupMember *p = static_cast<ThreadGroupMember *>( arg );
    p->c4_run_thread();
    return NULL;
}

//---------------------------------------------------------------------------//
// Kick off a collection of threads.
//---------------------------------------------------------------------------//

template<class T>
ThreadGroup<T>::ThreadGroup( int nthreads_, const T& model /*= T()*/ )
    : nthreads( nthreads_ ), v( nthreads )
{
    tcb = new ThreadControl( nthreads );

    for( int i=0; i < nthreads; i++ ) {
        v[i] = new T( model );
        ThreadGroupMember *p = v[i];
        p->thread = i;
        p->threads = nthreads;
        p->tcb = tcb;

        int rc = pthread_create( &p->tid, NULL, &thread_start_routine, p );
    }
}

//---------------------------------------------------------------------------//
// Reap the collection of threads, deallocate control structure for the
// group, etc.
//---------------------------------------------------------------------------//

template<class T>
ThreadGroup<T>::~ThreadGroup()
{
    using namespace std;

    cout << "Reaping threads.\n";

    for( int i=0; i < nthreads; i++ )
        if (pthread_join( v[i]->tid, NULL ))
            cout << "Failed reaping " << i << endl;
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of ThreadGroup.cc
//---------------------------------------------------------------------------//
