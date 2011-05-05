//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTermination_Detector.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "../ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "../Termination_Detector.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

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
        {
            ut.failure("Termination_Detection did NOT detect nontermination.");
        }
        else
        {
            ut.passes("Termination_Detection detected nontermination.");
        }
    }

    // Will hang if the unit fails.  Unfortunately, there's no other
    // portable way to test.
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
//                   end of tstTermination_Detector.cc
//---------------------------------------------------------------------------//
