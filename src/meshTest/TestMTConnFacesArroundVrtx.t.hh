//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTConnFacesArroundVrtx.t.hh
 * \author Randy M. Roberts
 * \date   Fri Sep 10 12:55:42 1999
 * \brief  Implementation file for the TestMTConnFacesArroundVrtx class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestMTConnFacesArroundVrtx.hh"
#include "StructuredMeshCFAVTest.hh"
#include "UnStructuredMeshCFAVTest.hh"

#include <iostream>
using std::endl;

namespace rtt_meshTest
{

template<class MTFactory>
void TestMTConnFacesArroundVrtx<MTFactory>::error(bool &passed,
					     const std::string &msg)
{
    if (!passed)
    {
	os_m << "TestMTConnFacesArroundVrtx failed: " << msg << endl;
	passed_m = false;
    }

    // reset the variable
    passed = true;
}

template<class MTFactory>
void TestMTConnFacesArroundVrtx<MTFactory>::run()
{
    passed_m = true;
    run(Structuring());
}

template<class MTFactory>
void TestMTConnFacesArroundVrtx<MTFactory>::run(Structured)
{
    os_m << "In Structured tests."
	 << endl;

    StructuredMeshCFAVTest<MTFactory> tester(meshFactory_m);
    tester.run();

    typedef typename StructuredMeshCFAVTest<MTFactory>::MsgList MsgList;

    // Go through the messages from the tester, and find print out the
    // failures.
    
    for (MsgList::const_iterator mlit = tester.msgList().begin();
	 mlit != tester.msgList().end();
	 mlit++)
    {
	bool passed;
	passed = mlit->first;
	if (passed)
	    os_m << mlit->second << endl;
	else
	    error(passed, mlit->second);
    }
}

template<class MTFactory>
void TestMTConnFacesArroundVrtx<MTFactory>::run(UnStructured)
{
    os_m << "In UnStructured tests."
	 << endl;

    UnStructuredMeshCFAVTest<MTFactory> tester(meshFactory_m);
    tester.run();

    typedef typename UnStructuredMeshCFAVTest<MTFactory>::MsgList MsgList;

    // Go through the messages from the tester, and find print out the
    // failures.
    
    for (MsgList::const_iterator mlit = tester.msgList().begin();
	 mlit != tester.msgList().end();
	 mlit++)
    {
	bool passed;
	passed = mlit->first;
	error(passed, mlit->second);
    }
}

} // end namespace rtt_meshTest


//---------------------------------------------------------------------------//
//                              end of TestMTConnFacesArroundVrtx.t.hh
//---------------------------------------------------------------------------//
