//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstParallelUnitTest.cc
 * \author Kent Budge
 * \date   Thu Jun 1 17:42:58 2006
 * \brief  Test the functionality of the class ParallelUnitTest
 * \note   © Copyright 2006 Los Alamos National Securities, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>

#include "../ParallelUnitTest.hh"
#include "ds++/Release.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstOne( UnitTest &unitTest )
{
    unitTest.passes("Looks like the PASSMSG macro is working.");
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for ParallelUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        ParallelUnitTest ut( argc, argv, release );
        tstOne(ut);

        ut.status();
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
//                   end of tstunit_test.cc
//---------------------------------------------------------------------------//
