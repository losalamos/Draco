//----------------------------------*-C++-*----------------------------------//
// tSyncAll.cc
// Maurice LeBrun
// Wed Jan 25 14:32:04 1995
//---------------------------------------------------------------------------//
// @> Test case for all of the block-sync or block-spinlock combinations.
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "c4/Sync.hh"
#include "c4/SpinLock.hh"
using namespace C4;

#include <stdlib.h>
#include <iostream.h>

int do_sync_test(int sync_type);

//---------------------------------------------------------------------------//
// Tests Sync, Spinlock, and all derived synchronization classes.  Specify on
// the command line the number of the test you want to run, or nothing if you
// want them all run.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    int sync_type = (argc > 1) ? atoi(argv[1]) : -1;
#ifndef __KCC
    ios::sync_with_stdio();
#endif

    C4_Init( argc, argv );

    if (sync_type >= 0)
	do_sync_test(sync_type);
    else {
	for (int i = 0; ; i++) {
	    if (do_sync_test(i))
		break;
	}
    }

    C4_Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
// Function that does all the actual work.
//---------------------------------------------------------------------------//

#define DOBLOCK(SyncType) \
    if (C4_node() == 0) \
	cout << "Testing: " << #SyncType << "\n\n" << flush; \
    gsync(); \
    cout << "Hello #1 from " << node() << endl << flush; \
    { \
	SyncType sync; \
	cout << "Hello #2 from " << node() << endl << flush; \
    } \
    cout << "Hello #3 from " << node() << endl << flush; \
    gsync(); \
    if (node() == 0) \
	cout << "\nEnd of test\n\n" << flush;

int do_sync_test(int sync_type)
{
    switch (sync_type) {
    case 0:
	DOBLOCK(HSync);
	break;
    case 1:
	DOBLOCK(TSync);
	break;
    case 2:
	DOBLOCK(HTSync);
	break;
    case 3:
	DOBLOCK(SpinLock);
	break;
    case 4:
	DOBLOCK(HSyncSpinLock);
	break;
    case 5:
	DOBLOCK(TSyncSpinLock);
	break;
    case 6:
	DOBLOCK(HTSyncSpinLock);
	break;
    default:
	return 1;
    }
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tSyncAll.cc
//---------------------------------------------------------------------------//
