//----------------------------------*-C++-*----------------------------------//
// tBandMat.cc
// Geoffrey M. Furnish
// Fri Mar 13 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "linalg/Banded_Matrix.cc"

#include "c4/global.hh"
#include "c4/Baton.hh"
#include "c4/SpinLock.hh"
using namespace C4;

#include <iostream>
using namespace std;

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

struct layout
{
    int nro, nrp, nrt;

    int row_offset() const { return nro; }
    int nrows_this_processor() const { return nrp; }
    int nrows_total() const { return nrt; }
//     bool verbose() const { return true; }
    bool verbose() const { return false; }
};

//---------------------------------------------------------------------------//
// Calculate the row decomposition for the problem.
//---------------------------------------------------------------------------//

void setup_layout( layout& l, int nrt )
{
    cout << "setting up layout on " << node() << endl;
    int nrp = nrt / nodes();
    if (node() < nrt - nrp*nodes())
	nrp++;
    int x = nrp;
    gsum(x);
    Insist( x == nrt, "Can't count!" );

    l.nrp = nrp;
    l.nrt = nrt;

    {
	Baton<int> b(0);
	l.nro = b;
	b += nrp;
    }

    cout << flush;
    {
	HTSyncSpinLock h;

	cout << node() << " l.nro=" << l.nro
	     << " l.nrp=" << l.nrp
	     << " l.nrt=" << l.nrt
	     << endl;
    }
}

void t1()
{
    int doff[3] = {-6, 0, 6};
    layout l;
    setup_layout( l, 8 );

    Banded_Matrix<int,3> m( l, doff );

// Build up a 3-diagonal identity matrix.  And you wonder why we need 3
// diagonals to test the identity multiplyer :-)...

    for( int i=0; i < l.nrp; i++ ) {
	m(0,i) = 0;
	m(1,i) = 1;
	m(2,i) = 0;
    }

    Mat1<int> x(l.nrp), b(l.nrp);

// Initialize x;
    for( int i=0; i < l.nrp; i++ ) {
	x(i) = l.nro + i;
	b(i) = 99;
    }

    m.multiply( x, b );

    cout << flush; gsync();

// Compare b with known correct answer.

    {
	HTSyncSpinLock h;

	for( int i=0; i < l.nrp; i++ )
	    cout << "nro+i=" << l.nro+i
		 << "  x(i)=" << x(i)
		 << "  b(i)=" << b(i) << endl;
    }

    int sf = 1;
    for( int i=0; i < l.nrp; i++ )
	if (x(i) != b(i)) sf = 0;
    gsum(sf);

    if (node() == 0)
	report( "t1", sf == nodes() );
}

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    if (node() == 0)
	cout << "Testing the Banded_Matrix clas.\n";

    try {
	t1();
    }
    catch( assertion& a )
    {
	cout << "Node " << node() << " caught assertion: " << a.what() << endl;
    }

    if (node() == 0)
	cout << "Testing complete.\n";

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tBandMat.cc
//---------------------------------------------------------------------------//
