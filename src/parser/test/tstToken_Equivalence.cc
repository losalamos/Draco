//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstToken_Equivalence.cc
 * \author Kelly Thompson
 * \date   Fri Jul 21 09:10:49 2006
 * \brief  Unit test for functions in Token_Equivalence.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include "ds++/ScalarUnitTest.hh"
#include "../Token_Equivalence.hh"
#include "ds++/Release.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstOne( UnitTest & ut )
{
    // create some data
    ostringstream data;
    data << "METIS decomposition specified\n"
         << "Dump cycle interval defaulting to 1.\n"
         << "Cycle       : 0\n"
         << "Time Step   : 1e-16 s.\n"
         << "Problem Time: 0 s.\n"
         << "error(0): 1     spr: 1\n"
         << "error(1): 0.00272636     spr: 0.00272636\n"
         << "error(2): 8.14886e-06     spr: 0.00298892\n"
         << "pid[0] done error(2): 8.14886e-06  spr: 0.00298892\n"
         << "User Cpu time this time step: 5.17 \n" << endl;

    // create a string token stream from the data

    String_Token_Stream tokens( data.str() );

    // Test Token_Equivalence functions:

    // look for an int associated with a keyword
    check_token_keyword_value( tokens, "Cycle", 0, ut );

    // look for a double associated with a keyword (4th occurance).
    check_token_keyword_value( tokens, "spr", 0.00298892, ut, 4 );

    // look for a keyword.
    check_token_keyword( tokens, "User Cpu time this time step", ut );
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        ScalarUnitTest ut( argc, argv, release );
        tstOne(ut);
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
//                        end of tstToken_Equivalence.cc
//---------------------------------------------------------------------------//
