//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTermination_Detector.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \brief  
 * \note   Copyright (C) 2004-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Termination_Detector.hh"
#include "../ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <iostream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

#define PASSMSG(A) ut.passes(A)
#define FAILMSG(A) ut.failure(A)
#define ITFAILS    ut.failure( __LINE__ )

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstTermDet( UnitTest & ut )
{
    
    Termination_Detector td(1);

    td.init();

    td.update_receive_count(0);
    td.update_send_count(1);
    td.update_work_count(2);

    for (unsigned c=0; c<5; ++c)
    {
        if( td.is_terminated() )
            FAILMSG("Termination_Detection did NOT detect nontermination.");
        else
            PASSMSG("Termination_Detection detected nontermination.");
    }

    // Will hang if the unit fails.  Unfortunately, there's no other portable
    // way to test.
    td.update_receive_count(1);

    while (!td.is_terminated()) {/* do nothing */};
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    try
    {
       ParallelUnitTest ut( argc, argv, release );
       tstTermDet(ut);
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown" << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstTermination_Detector.cc
//---------------------------------------------------------------------------//
