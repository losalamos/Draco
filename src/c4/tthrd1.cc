//----------------------------------*-C++-*----------------------------------//
// tthrd1.cc
// Geoffrey Furnish
// Mon Jul 13 22:48:24 1998
//---------------------------------------------------------------------------//
// @> Test ThreadGroup and ThreadGroupMember classes.
//---------------------------------------------------------------------------//

#include "c4/ThreadGroup.cc"
#include "c4/ThreadGroupMember.hh"

#include <iostream>

void z( int& x )
{
    x++;
}

class X : public C4::ThreadGroupMember
{
  public:
    void c4_run_thread()
    {
        int x = 1000000 * (rand() & 7);
        int y = 0;

        for( int i=0; i < x; i++ )
            z( y );

        printf( "Hello #1 from thread member %d of %d.\n", thread, threads );

        gsync();

        int s = thread;
        gsum(s);
        printf( "Thread %d, sum of thread id's is %d\n", thread, s );

        s = thread+101;
        gmin(s);
        printf( "Thread %d, min of operands is %d\n", thread, s );

        x = 1000000 * (rand() & 7);

        for( int i=0; i < x; i++ )
            z( y );

        printf( "Hello #2 from thread member %d of %d.\n", thread, threads );
    }
};

int main( int argc, char *argv[] )
{
    using namespace  std;
    using namespace  C4;

    cout << "Initiating thread test program.\n";

    {
        ThreadGroup<X> tgx( 5 );
        cout << "Thread group is active.\n";
    }

    fflush( stdout );

    cout << "Thread group has gone out of scope.\n";
}

//---------------------------------------------------------------------------//
//                              end of tthrd1.cc
//---------------------------------------------------------------------------//
