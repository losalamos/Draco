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
void TestMTConnFacesArroundVrtx<MTFactory>::run()
{
    setPassed(true);
    run(Structuring());
}

template<class MTFactory>
void TestMTConnFacesArroundVrtx<MTFactory>::run(Structured)
{
    os() << "In Structured tests."
	 << endl;

    StructuredMeshCFAVTest<MTFactory> tester(*this, meshFactory_m);
    tester.run();
}

template<class MTFactory>
void TestMTConnFacesArroundVrtx<MTFactory>::run(UnStructured)
{
    os() << "In UnStructured tests."
	 << endl;

    UnStructuredMeshCFAVTest<MTFactory> tester(*this, meshFactory_m);
    tester.run();
}

} // end namespace rtt_meshTest


//---------------------------------------------------------------------------//
//                              end of TestMTConnFacesArroundVrtx.t.hh
//---------------------------------------------------------------------------//
