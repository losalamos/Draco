//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstConstants.cc
 * \author Todd J. Urbatsch
 * \date   Fri Apr 11 09:10:31 2003
 * \brief  Constants.hh tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "../Constants.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <cmath>

using rtt_dsxx::soft_equiv;

using std::cout;
using std::endl;
using std::string;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// test the constants
//---------------------------------------------------------------------------//
/*!  
 * \brief Tests the constants in mc/Constants.hh, some of which are derived.
 */
void test_constants()
{
    using rtt_mc::global::pi;
    using rtt_mc::global::huge;
    using rtt_mc::global::huge_int;
    using rtt_mc::global::epsilon;
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using rtt_mc::global::k;
    using rtt_mc::global::h;

    if (!soft_equiv(pi, 3.1415926, 1.0e-7))   ITFAILS;
    if (huge < 1.0e7)                         ITFAILS;
    if (huge_int < 10000000)                  ITFAILS;
    if (epsilon > 0.01)                       ITFAILS;
    if (!soft_equiv(a, 0.01372, 1.0e-4))      ITFAILS;
    if (!soft_equiv(c, 299.792, 1.0e-5))      ITFAILS;
    if (!soft_equiv(k, 8.617308e-8, 1.0e-6))  ITFAILS;
    if (!soft_equiv(h, 4.135706e-18, 1.0e-6)) ITFAILS;

}
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
        C4::Finalize();
        return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

        // test the constants
	test_constants();

    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstConstants, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstConstants Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstConstants on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                   end of test/tstConstants.t.hh
//---------------------------------------------------------------------------//
