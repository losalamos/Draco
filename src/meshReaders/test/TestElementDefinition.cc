//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestElementDefinition.cc
 * \author John McGhee
 * \date   Fri Mar  3 08:41:46 2000
 * \brief  Implements the Element_Definition unit tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestElementDefinition.hh"
#include "../Element_Definition.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>

using std::cout;
using std::endl;


namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_meshReaders_test::TestElementDefinition;
    
    return SP<TestApp>(new TestElementDefinition(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_meshReaders_test
{

using std::string;
using rtt_meshReaders::Element_Definition;

TestElementDefinition::TestElementDefinition(int argc, char *argv[],
					       std::ostream& os_in)

    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestElementDefinition" << endl;
}

string TestElementDefinition::version() const
{
    return rtt_meshReaders::release();
}

string TestElementDefinition::runTest()
{
    fail() << "This test sucks!";

    if (passed())
    {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestElementDefinition.cc
//---------------------------------------------------------------------------//
