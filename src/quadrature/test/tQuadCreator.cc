//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tQuadCreator.cc
 * \author Kelly Thompson
 * \date   Fri Aug 18 12:46:41 2006
 * \brief  Unit test for quadcreator class.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//


#include <string>
#include <iostream>
#include <sstream>

#include "ds++/SP.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/String_Token_Stream.hh"

#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "ds++/Release.hh"

using namespace rtt_dsxx;
using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstTokenCtor1( rtt_dsxx::ScalarUnitTest & ut )
{
    using namespace std;
    using namespace rtt_quadrature;

    cout << "\n>>> tstTokenCtor1 <<<\n" << endl;
    
    try
    {
        ostringstream contents;

        contents << "type = level symmetric\n"
                 << "order = 2\n"
                 << "end\n"
                 << endl;
        
        rtt_parser::String_Token_Stream tokens( contents.str() );
        
        // Create a Quadrature from this set of tokens.
        
        rtt_dsxx::SP< Quadrature const > spQuad = QuadCreator().quadCreate( tokens );
        
        if( spQuad->name() == "2D Level Symmetric" )
            ut.passes("Found correct quadrature name.");
        else
            ut.failure("Did not find expected quadrature name.");
        
        if( spQuad->getSnOrder() == 2 )
            ut.passes("Found correct SnOrder.");
        else
            ut.failure("Did not find expected SnOrder.");

        // test write/read

        string const text = spQuad->as_text("\n");

        tokens = String_Token_Stream(text);

        SP<Quadrature const> spTextWuad = QuadCreator().quadCreate(tokens);

        if (tokens.error_count()==0)
        {
            ut.passes("write\reawd is correct");
        }
        else
        {
            ut.failure("write\read is NOT correct");
        }
    }
    catch(...)
    {
        ut.failure("Encountered a C++ Exception.");
        throw;
    }
    return;
}

//---------------------------------------------------------------------------//

void tstTokenCtor2( rtt_dsxx::ScalarUnitTest & ut )
{
    using namespace std;
    using namespace rtt_quadrature;

    cout << "\n>>> tstTokenCtor2 <<<\n" << endl;

    try
    {
        ostringstream contents;

        // All caps set to lc.
        // order must be even.
        
        contents << "sQuaRe CL\n"
                 << "order 16\n"
                 << "end\n"
                 << endl;

        rtt_parser::String_Token_Stream tokens( contents.str() );
        
        // Create a Quadrature from this set of tokens.
        
        rtt_dsxx::SP< Quadrature const > spQuad = QuadCreator().quadCreate( tokens );

        if( spQuad->name() == "2D Square Chebyshev Legendre" )
            ut.passes("Found correct quadrature name.");
        else
            ut.failure("Did not find expected quadrature name.");
        
        if( spQuad->getSnOrder() == 16 )
            ut.passes("Found correct SnOrder.");
        else
            ut.failure("Did not find expected SnOrder.");
        
    }
    catch(...)
    {
        ut.failure("Encountered a C++ Exception.");
        throw;
    }
    return;
}
//---------------------------------------------------------------------------//

void tstTokenCtor3( rtt_dsxx::ScalarUnitTest & ut )
{
    using namespace std;
    using namespace rtt_quadrature;

    cout << "\n>>> tstTokenCtor3 <<<\n" << endl;
    
    try
    {
        ostringstream contents;

        contents << "type = gauss legendre\n"
                 << "order = 128\n"
                 << "interpolation algorithm = Galerkin\n"
                 << "end\n"
                 << endl;

        rtt_parser::String_Token_Stream tokens( contents.str() );
        
        // Create a Quadrature from this set of tokens.
        
        rtt_dsxx::SP< Quadrature const > spQuad = QuadCreator().quadCreate( tokens );

        if( spQuad->name() == "1D Gauss Legendre" )
            ut.passes("Found correct quadrature name.");
        else
            ut.failure("Did not find expected quadrature name.");
        
        if( spQuad->getSnOrder() == 128 )
            ut.passes("Found correct SnOrder.");
        else
            ut.failure("Did not find expected SnOrder.");
        
    }
    catch(...)
    {
        ut.failure("Encountered a C++ Exception.");
        throw;
    }
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using namespace std;

    try
    {
        rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
        tstTokenCtor1( ut );
        tstTokenCtor2( ut );
        tstTokenCtor3( ut );
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
//                        end of tQuadrature.cc
//---------------------------------------------------------------------------//
