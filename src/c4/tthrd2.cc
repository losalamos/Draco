//----------------------------------*-C++-*----------------------------------//
// tthrd2.cc
// Geoffrey Furnish
// Mon Aug 4 14:48:24 1998
//---------------------------------------------------------------------------//
// @> Test ThreadGroup and ThreadGroupMember classes.
//---------------------------------------------------------------------------//

#include "c4/ThreadGroup.cc"
#include "c4/ThreadGroupMember.hh"
#include "c4/global.hh"

#include <iostream>
#include <vector>
using namespace std;

int nthreads = 2;

class X : public C4::ThreadGroupMember
{
  public:
    void c4_run_thread()
    {
        vector<int> v;
        v.push_back( 1 );
        v.push_back( 10 );
        v.push_back( 100 );
        v.push_back( 1000 );
        v.push_back( 10000 );
        v.push_back( 100000 );
        v.push_back( 1000000 );

        int niters = 25;

        if (thread == 0)
            printf( "Timing reductions with %d threads.\n", threads );

        double start = C4::Wtime();

        for( int i=0; i < v.size(); i++ )
        {
            double s = t1( v[i], niters );
        //             double avg = s;
        //             C4::gsum(avg);
        //             avg /= nn;

            if (thread == 0)
                printf( "Elements: %7d, microseconds/pass  thread 0: %10.2lf\n",
                        v[i], s / 1.0e-6 );
        }

        double end = C4::Wtime();

        if (thread == 0)
            printf( "Total time to run was %lf seconds.\n", end-start );

    }

    double t1( int nitems, int niters )
    {
    // We make one extra buffer so we can do one untimed reduction.

        int **ppx = new int*[ niters+1 ];
        for( int i=0; i < niters+1; i++ )
        {
            ppx[i] = new int[ nitems ];

            for( int j=0; j < nitems; j++ )
                ppx[i][j] = j;
        }

    // Prime the pump by doing one reduction outside the timing loop.

        gsum( ppx[niters], nitems );

    // Now time a bunch of reduction operations.

        double t1 = C4::Wtime();

        for( int i=0; i < niters; i++ )
            gsum( ppx[i], nitems );

        double t2 = C4::Wtime();

    // Now check the results.

        int bogus=0;

        for( int i=0; i < niters+1; i++ )
            for( int j=0; j < nitems; j++ )
                if (ppx[i][j] != threads*j)
                    bogus++;

        if (bogus) {
            lock();

            printf( "Thread %d detected %d bogus results.\n", thread, bogus );

            for( int j=0; j < nitems; j++ )
            {
                printf( "j=%4d ", j );
                for( int i=0; i < 1; i++ )
                {
                    printf( " %7d", ppx[i][j] );
                }
                printf( "\n" );
            }

            unlock();
        }

    // Now clean up.

        for( int i=0; i < niters+1; i++ )
            delete[] ppx[i];
        delete[] ppx;

        return (t2 - t1) / niters;
    }
};

int main( int argc, char *argv[] )
{
    using namespace  std;
    using namespace  C4;

    cout << "Testing (timing) array thread reductions.\n";

    if (argc == 2)
        nthreads = atoi( argv[1] );

    {
        ThreadGroup<X> tgx( nthreads );
    }

    fflush( stdout );
}

//---------------------------------------------------------------------------//
//                              end of tthrd2.cc
//---------------------------------------------------------------------------//
