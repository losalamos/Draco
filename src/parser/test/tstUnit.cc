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
#include "../Release.hh"
#include "../Unit.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"

using namespace std;

using namespace rtt_parser;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void unit_test()
{
    Unit tstC  = { 0, 0,  1, 1, 0, 0, 0, 0, 0,    1};
    if (C!=tstC)
    {
	FAILMSG("unit C does NOT have expected dimensions");
    }
    else
    {
	PASSMSG("unit C has expected dimensions");
    }

    Unit tstHz = { 0, 0, -1, 0, 0, 0, 0, 0, 0,    1};
    if (Hz!=tstHz)
    {
	FAILMSG("unit Hz does NOT have expected dimensions");
    }
    else
    {
	PASSMSG("unit Hz has expected dimensions");
    }

    Unit tstN  = { 1, 1, -2, 0, 0, 0, 0, 0, 0,    1}; 
    if (N!=tstN)
    {
	FAILMSG("unit N does NOT have expected dimensions");
    }
    else
    {
	PASSMSG("unit N has expected dimensions");
    }

    Unit tstJ  = { 2, 1, -2, 0, 0, 0, 0, 0, 0,    1};
    if (J!=tstJ)
    {
	FAILMSG("unit J does NOT have expected dimensions");
    }
    else
    {
	PASSMSG("unit J has expected dimensions");
    }
    {
        ostringstream buffer;
        buffer << tstJ;
        if (buffer.str() == "1 m^2-kg-s^-2")
        {
            PASSMSG("correct text representation of J");
        }
        else
        {
            FAILMSG("NOT correct text representation of J");
        }
    }

    Unit tstinch = {1, 0, 0, 0, 0, 0, 0, 0, 0,    0.0254};
    if (inch == tstinch)
    {
	PASSMSG("unit inch has expected dimensions");
    }
    else
    {
	FAILMSG("unit inch does NOT have expected dimensions");
    }

    // Test divisor
    {
        Unit five_cm = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.05};
        Unit one_cm  = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.01};
        if( five_cm/5.0 == one_cm )
        {
            PASSMSG("Units: Division by a constant works.");
        }
        else
        {
            ostringstream buffer;
            buffer << "Units: Division by a constant fails.\n\t"
                   << "five_cm/5.0 = " << five_cm/5.0 << "\n\t"
                   << "one_cm = " << one_cm << "\n\t"
                   << "test five_cm/5.0 == one_cm ?" << std::endl;
            FAILMSG( buffer.str() );
        }
    }
    
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_parser::release() 
		 << endl;
            rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	if ( rtt_c4::nodes() == 1 ) unit_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstUnit, " << ass.what()
	     << endl;
        rtt_c4::finalize();
	return 1;
    }

    // status of test
    cout <<     "\n*********************************************\n";
    if( rtt_parser_test::passed )
        cout << "**** tstUnit Test: PASSED\n";
    cout <<     "*********************************************\n" << endl;

    rtt_c4::global_barrier();
    cout << "Done testing tstUnit." << endl;
    rtt_c4::finalize();
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstUnit.cc
//---------------------------------------------------------------------------//
