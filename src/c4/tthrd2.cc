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

int nthreads = 2;
int nitems = 100;
int niters = 100;

void z( int& x )
{
    x++;
}

class X : public C4::ThreadGroupMember
{
  public:
    void c4_run_thread()
    {
        int *px = new int[ nitems ];

        for( int i=0; i < nitems; i++ )
            px[i] = i;

        double t1 = C4::Wtime();

        for( int j=0; j < niters; j++ )
            gsum( px, nitems );

        double t2 = C4::Wtime();

        printf( "Thread %d took %f microseconds/reduction of %d elements\n",
                thread, (t2 - t1) / niters / 1.0e-6, nitems );
    }
};

int main( int argc, char *argv[] )
{
    using namespace  std;
    using namespace  C4;

    cout << "Testing (timing) array thread reductions.\n";

    if (argc == 4)
    {
        nthreads = atoi( argv[1] );
        nitems = atoi( argv[2] );
        niters = atoi( argv[3] );
    }

    {
        ThreadGroup<X> tgx( nthreads );
    }

    fflush( stdout );
}

//---------------------------------------------------------------------------//
//                              end of tthrd2.cc
//---------------------------------------------------------------------------//
