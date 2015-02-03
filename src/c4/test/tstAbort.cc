//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4_test/test/tstAbort.cc
 * \author Thomas M. Evans
 * \date   Thu Jun  2 09:28:02 2005
 * \brief  C4 Abort test.
 * \note   Copyright (C) 2005-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ParallelUnitTest.hh"
#include "ds++/Release.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void abort_test(rtt_dsxx::UnitTest & ut)
{
    using namespace std;
    
    cout << "Entering abort_test on proc " << rtt_c4::node() << endl;
    
    rtt_c4::global_barrier();

    // only abort from processor 0 for nice output
    if (rtt_c4::node() == 0)
    {
        cout << "Aborting from processor 0" << endl;
        rtt_c4::abort();
        FAILMSG("Should never get here.");
    }
    else
    {
        PASSMSG("Only abort from Processor 0.");
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);

    // runtest command option tag
    bool runtest = false;
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[ arg ])=="--runtest") runtest = true;

    try
    { 
        // run test here so we get a pass message; this should simply abort the
        // program at this point;  only run the test if --runtest is given so we
        // don't hose the automated testing
        if (runtest)
        {
            abort_test(ut);
            std::cout<<"Do we get here?"<<std::endl;
        }
        else
            PASSMSG("Test allways passes when --runtest is not provided.");
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstAbort.cc
//---------------------------------------------------------------------------//
