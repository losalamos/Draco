//----------------------------------*-C++-*----------------------------------//
// tsr3.cc
// Geoffrey M. Furnish
// Mon Jun 22 08:11:39 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../global.hh"

#include "ds++/Mat.hh"

#include <iostream>
using namespace std;

void t1( int len, int msgs )
{
    using namespace dsxx;

    Mat1<int> buf1( len ), buf2( len );

    int node = C4::node();

    if (node == 0)
    {
        cout << "Sending " << msgs
             << " messages of length " << len*sizeof(int) << endl;

        for( int i=0; i < len; i++ )
            buf1(i) = i;

        double start = C4::Wtime();

        for( int i=0; i < msgs; i++ )
        {
            C4::Send( &buf1(0), len, 1 );
            C4::Recv( &buf2(0), len, 1 );
        }

        double end = C4::Wtime();

        double seconds = end - start;
        double microsecs = seconds / 1.0e-6;

        cout << "It took " << end - start << " seconds.\n";
        cout << "Round trip time per message is "
             << microsecs / msgs << " microseconds/message.\n";
    }
    if (node == 1)
    {
        for( int i=0; i < len; i++ )
            buf2(i) = i;

        for( int i=0; i < msgs; i++ )
        {
            C4::Recv( &buf1(0), len, 0 );
            C4::Send( &buf2(0), len, 0 );
        }
    }
}

int size = 1000;
int cnt = 500;

void process_cli( int argc, char  *argv[] )
{
    argc--, argv++;             // skip program name

    while( argc )
    {
        if (!strcmp(*argv,"-s"))
        {
            size = atoi( argv[1] );
            argc -= 2, argv += 2;
            continue;
        }

        if (!strcmp(*argv,"-c"))
        {
            cnt = atoi( argv[1] );
            argc -= 2, argv += 2;
            continue;
        }

        cout << "unrecognized argument " << *argv << endl;
        argc--, argv++;
    }
}

int main( int argc, char  *argv[] )
{
    process_cli( argc, argv );

    C4::Init( argc, argv );

    C4::gsync();

    t1( size, cnt );

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tsr3.cc
//---------------------------------------------------------------------------//
