//----------------------------------*-C++-*----------------------------------//
// main.cc
// Geoffrey Furnish
// Thu Jun 26 21:08:23 1997
//---------------------------------------------------------------------------//
// @> Draco bootstrap vehicle.
//---------------------------------------------------------------------------//

#include <iostream.h>

#include "c4/global.hh"
#include "c4/SpinLock.hh"

#include "Pooma.h"

#include "draco.hh"

int main( int argc, char *argv[] )
{
    Pooma pooma( argc, argv );

    cout << "Draco is running.\n";

    {
	HTSyncSpinLock h;
        cout << "nodes = " << pooma.getNodes() << endl;
	cout << "Hello from node " << C4_node() << endl;
    }

    {
	draco d;

	cout << "ready to fetch a k_mesh from draco." << endl;

    }

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of main.cc
//---------------------------------------------------------------------------//
