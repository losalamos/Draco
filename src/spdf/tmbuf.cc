//----------------------------------*-C++-*----------------------------------//
// tmbuf.cc
// Geoffrey Furnish
// Fri Aug 11 14:28:09 1995
//---------------------------------------------------------------------------//
// @> Test the mbuf_spdf_stream.
//---------------------------------------------------------------------------//

#include <iostream.h>

#include "Assert.hh"

#include "spdf/spdf_stream.hh"

int main()
{
    int i, len;
    mbuf_spdf_stream ms;

    float *fi = new float[10];
    float *fo = new float[10];

    for( i=0; i < 10; i++ )
	fi[i] = i*i;

    ms.wrp( fi, 10 );
    ms.reset();
    ms.rdp( fo, len, 10 );

    cout << "Got back " << len << ", expected 10.\n";

    for( i=0; i < len; i++ )
	cout << "in[" << i << "]=" << fi[i] << " out=" << fo[i] << endl;

// Now let's test the double write packed stuff.

    double *di = new double[7];
    double *od = new double[7];	// "do" is reserved word, grrrr.

    for( i=0; i < 7; i++ )
	di[i] = 500000. + i*i*i*i;

    ms.reset();
    cout << "writing packed double" << endl;
    ms.wrp( di, 7 );
    ms.reset();
    cout << "reading packed double" << endl;
    ms.rdp( od, len, 7 );

    cout << "Reading double packed, got " << len
	 << " expected 7." << endl;

    for( i=0; i < 7; i++ )
	cout << "in[" << i << "]= " << di[i]
	     << " out=" << od[i] << endl;
}

//---------------------------------------------------------------------------//
//                              end of tmbuf.cc
//---------------------------------------------------------------------------//
