#include "../global.hh"

#include <iostream>
#include <vector>
using namespace std;

//---------------------------------------------------------------------------//
// This function actually performs the global reduction and produces the
// average time per reduction on each node.
//---------------------------------------------------------------------------//

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

    C4::gsum( ppx[niters], nitems );

// Now time a bunch of reduction operations.

    double t1 = C4::Wtime();

    for( int i=0; i < niters; i++ )
        C4::gsum( ppx[i], nitems );

    double t2 = C4::Wtime();

// Now clean up.

    for( int i=0; i < niters+1; i++ )
        delete[] ppx[i];
    delete[] ppx;

    return (t2 - t1) / niters;
}

//---------------------------------------------------------------------------//
// The purpose of this program is to collect timing data on global reduction
// routines in MPI and SHMEM.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
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

    C4::Init( argc, argv );

    int nn = C4::nodes();

    if (C4::node() == 0)
        printf( "Timing reductions on %d nodes.\n", nn );

    double start = C4::Wtime();

    for( int i=0; i < v.size(); i++ )
    {
        double s = t1( v[i], niters );
        double avg = s;
        C4::gsum(avg);
        avg /= nn;
        if (C4::node() == 0)
            printf( "Elements: %7d, microseconds/pass  node 0: %lf  avg: %lf\n",
                    v[i], s / 1.0e-6, avg / 1.0e-6 );
    }

    double end = C4::Wtime();

    if (C4::node() == 0)
        printf( "Total time to run was %lf seconds.\n", end-start );

    C4::Finalize();
}
