//----------------------------------*-C++-*----------------------------------//
// ttest.cc
// Geoffrey M. Furnish
// Mon Jul 13 11:46:50 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../global.hh"
#include "../C4_Req.hh"

#include <iostream>

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    C4::gsync();

    if (C4::node() == 0)
    {
        double beg = C4::Wtime();

        while( C4::Wtime() - beg < 1.0 )
            ;                   // Kill some time.

        C4::Send( 7, 1, 42 );       // Send a 7 to node 1.
    }
    if (C4::node() == 1)
    {
        int x;
        C4::C4_Req r = C4::RecvAsync( &x, 1, 0, 42 );

        int cnt = 0;

        while( !r.complete() ) cnt++;

        std::cout << cnt
                  << " tests were conducted while waiting for completion.\n";
    }

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of ttest.cc
//---------------------------------------------------------------------------//
