//----------------------------------*-C++-*----------------------------------//
// tsr2.cc
// Geoffrey M. Furnish
// Fri Mar  6 11:36:43 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../global.hh"
using namespace C4;

#include <pthread.h>

#include <iostream>
using namespace std;
#include <algorithm>

#include <unistd.h>
#include <stdlib.h>

//---------------------------------------------------------------------------//
// Indicate the success or failure of a test.
//---------------------------------------------------------------------------//

void report( const char *test, bool sf )
{
    if (sf)
	cout << test << ": passed\n";
    else
	cout << test << ": failed\n";
}

//---------------------------------------------------------------------------//
// Simple test that we can find out our nodal identity and perform a barier
// synchronization.
//---------------------------------------------------------------------------//

void t1()
{
    printf( "Hello from node %d of %d\n", node(), nodes() );
    gsync();
}

//---------------------------------------------------------------------------//
// Simple test of blocking send.
//---------------------------------------------------------------------------//

void t2()
{
    int x;
    if (node() == 0) {
	x = 44;
	C4::Send( x, 1 );
    }
    if (node() == 1) {
	C4::Recv( x, 0 );
	report( "t2", x == 44 );
    }
}

//---------------------------------------------------------------------------//
// Function to be run in a thread, which does a blocking send.
//---------------------------------------------------------------------------//

void *t3_thread_blocking_send( void *arg )
{
    int x = *(int *) arg;

    pthread_detach( pthread_self() );

    sleep( rand() & 3 );

    C4::Send( x, 1 );

    return NULL;
}

//---------------------------------------------------------------------------//
// Test use of C4::Send inside of threads.
//---------------------------------------------------------------------------//

void t3()
{
    int nthreads = 200;

    if (node() == 0)
    {
	pthread_t thread;
	int *pids = new int[ nthreads ];

    // Spawn oodles of threads which xmit data to node 1.
	for( int i=0; i < nthreads; i++ ) {
	    pids[i] = i;
	    pthread_create( &thread, NULL,
			    t3_thread_blocking_send, &pids[i] );
	}
    }
    if (node() == 1)
    {
	int *pvals = new int[ nthreads ];

    // Receive all the messages sent by the threads on node 0.
	for( int i=0; i < nthreads; i++ )
	    C4::Recv( pvals[i], 0 );

    // Now validate them all.
	std::sort( pvals, pvals+nthreads );

	bool sf = true;
	for( int i=0; i < nthreads; i++ )
	    if (pvals[i] != i) sf = false;

    // Indicate status.
	report( "t3", sf );
    }
}

//---------------------------------------------------------------------------//
// The main function just runs the tests above.
//---------------------------------------------------------------------------//

void main( int argc, char *argv[] )
{
    Init( argc, argv );

    t1();
    t2();
    t3();

    Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tsr2.cc
//---------------------------------------------------------------------------//
