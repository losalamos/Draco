//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstUnit.cc
 * \author Kent G. Budge
 * \date   Feb 18 2003
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <sstream>
#include "parser_test.hh"
#include "ds++/Release.hh"
#include "../Unit.hh"
#include "ds++/Assert.hh"
#include "ds++/ScalarUnitTest.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

using namespace std;

using namespace rtt_parser;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void unit_test(UnitTest &ut)
{
    Unit tstC  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
    if (C!=tstC)
    {
	ut.failure("unit C does NOT have expected dimensions");
    }
    else
    {
	ut.passes("unit C has expected dimensions");
    }

    Unit tstHz = { 0, 0, -1, 0, 0, 0, 0, 0, 0,    1};
    if (Hz!=tstHz)
    {
	ut.failure("unit Hz does NOT have expected dimensions");
    }
    else
    {
	ut.passes("unit Hz has expected dimensions");
    }

    Unit tstN  = { 1, 1, -2, 0, 0, 0, 0, 0, 0,    1}; 
    if (N!=tstN)
    {
	ut.failure("unit N does NOT have expected dimensions");
    }
    else
    {
	ut.passes("unit N has expected dimensions");
    }

    Unit tstJ  = { 2, 1, -2, 0, 0, 0, 0, 0, 0,    1};
    if (J!=tstJ)
    {
	ut.failure("unit J does NOT have expected dimensions");
    }
    else
    {
	ut.passes("unit J has expected dimensions");
    }
    {
        ostringstream buffer;
        buffer << tstJ;
        if (buffer.str() == "1 m^2-kg-s^-2")
        {
            ut.passes("correct text representation of J");
        }
        else
        {
            ut.failure("NOT correct text representation of J");
        }
    }

    Unit tstinch = {1, 0, 0, 0, 0, 0, 0, 0, 0,    0.0254};
    if (inch == tstinch)
    {
	ut.passes("unit inch has expected dimensions");
    }
    else
    {
	ut.failure("unit inch does NOT have expected dimensions");
    }

    // Test divisor
    {
        Unit five_cm = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.05};
        Unit one_cm  = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.01};
        if( five_cm/5.0 == one_cm )
        {
            ut.passes("Units: Division by a constant works.");
        }
        else
        {
            ostringstream buffer;
            buffer << "Units: Division by a constant fails.\n\t"
                   << "five_cm/5.0 = " << five_cm/5.0 << "\n\t"
                   << "one_cm = " << one_cm << "\n\t"
                   << "test five_cm/5.0 == one_cm ?" << std::endl;
            ut.failure( buffer.str() );
        }
    }

    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0.1, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect length difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect length difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0.1,  1, 1, 0, 0, 0, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect mass difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect mass difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1.1, 1, 0, 0, 0, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect time difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect time difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1.1, 0, 0, 0, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect current difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect current difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1, 0.1, 0, 0, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect temperature difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect temperature difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1, 0, 0.1, 0, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect mole difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect mole difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1, 0, 0, 0.1, 0, 0,    1};
        if (t1==t2) ut.failure("failed to detect luminence difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect luminence difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1, 0, 0, 0, 0.1, 0,    1};
        if (t1==t2) ut.failure("failed to detect rad difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect rad difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1, 0, 0, 0, 0, 0.1,    1};
        if (t1==t2) ut.failure("failed to detect solid angle difference");
        if (is_compatible(t1, t2))
            ut.failure("failed to detect solid angle difference");
    }
    {
        Unit t1  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
        Unit t2  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1.1};
        if (t1==t2) ut.failure("failed to detect conversion factor difference");
        if (!is_compatible(t1, t2))
            ut.failure("failed to ignore conversion factor difference");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
	// >>> UNIT TESTS
        ScalarUnitTest ut( argc, argv, release );
        unit_test(ut);
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
//                        end of tstUnit.cc
//---------------------------------------------------------------------------//
