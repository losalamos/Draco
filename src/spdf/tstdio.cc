//----------------------------------*-C++-*----------------------------------//
// tstdio.cc
// Geoffrey Furnish
// Fri Aug 11 12:48:53 1995
//---------------------------------------------------------------------------//
// @> Test the stdio_spdf_stream;
//---------------------------------------------------------------------------//

#include <iostream.h>

#include "Assert.hh"

#include "spdf/spdf_stream.hh"

void arch();

void twrite();
void tread();

int main( int argc, char *argv[] )
{
    arch();

    Assert( argc == 2 );
    Assert( *argv[1] == 'r' || *argv[1] == 'w' || *argv[1] == 'd' );

    if (*argv[1] == 'w')
	twrite();
    else
	tread();

    return 0;
}

int iout[10] = { 1, 2, 3, 4, 5, 6, 7, 256, 65536, 1000000000 };
int iin[10];

float fout[5] = { 7., 7.e7, 7.e-7, -7.e7, -7.e-7 };
float fin[5];

float foutb[5] = { .15, 33.3, 48.8, 77.777, 99.9998789 };
float finb[5];

void arch()
{
    cout << "sizeof(int) = "    << sizeof(int)    << endl;
    cout << "sizeof(long) = "   << sizeof(long)   << endl;
    cout << "sizeof(long long) = "   << sizeof(long long)   << endl;
    cout << "sizeof(float) = "  << sizeof(float)  << endl;
    cout << "sizeof(double) = " << sizeof(double) << endl;

    cout << "Specializations:\n";
    cout << "spdf_traits<double>::shift_factor() = "
	 << spdf_traits<double>::shift_factor() << endl;
    cout << "sizeof(spdf_traits<double>::bitrep) = "
	 << sizeof(spdf_traits<double>::bitrep) << endl;
}

void twrite()
{
    stdio_spdf_stream s( "spdf.dat", "w" );

    s.wr( 7 );
    s.wr( 7.f );
    s.wr( 7. );

    s.wr( iout, 10 );
    s.wr( fout, 5 );

    s.wrp( foutb, 5 );
}

void tread()
{
    int i, ilen, j, k, flen;
    float f;
    double d;

    stdio_spdf_stream s( "spdf.dat", "r" );

    s.rd( i );
    s.rd( f );
    s.rd( d );

    cout << "i=" << i << " f=" << f << " d=" << d << endl;

    s.rd( iin, ilen, 10 );
    s.rd( fin, flen, 5 );

    Assert( ilen == 10 );
    for( j=0; j < 10; j++ )
	cout << "iin[" << j << "]=" << iin[j] << " "
	     << (iin[j] == iout[j] ? "ok" : "WRONG!") << endl;

    Assert( flen == 5 );
    for( k=0; k < 5; k++ )
	cout << "fin[" << k << "]=" << fin[k] << " "
	     << (fin[k] == fout[k] ? "ok" : "WRONG!") << endl;

    s.rdp( finb, flen, 5 );
    Assert( flen == 5 );
    for( k=0; k < 5; k++ )
	cout << "finb[" << k << "]=" << finb[k] << " "
	     << ( fabs(finb[k] - foutb[k])/foutb[k] < .01 ? "ok" : "WRONG!" )
	     << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstdio.cc
//---------------------------------------------------------------------------//
