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

struct layout
{
    int nro, nrp, nrt;

    int row_offset() const { return nro; }
    int nrows_this_processor() const { return nrp; }
    int nrows_total() const { return nrt; }
    bool verbose() const { return true; }
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
}

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    if (node() == 0)
	cout << "Testing the Banded_Matrix clas.\n";

    t1();

    if (node() == 0)
	cout << "Testing complete.\n";

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tBandMat.cc
//---------------------------------------------------------------------------//
