//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestHexFormat.cc
 * \author John McGhee
 * \date   Thu Mar  9 08:54:59 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestHexFormat.hh"
#include "../Hex_Format.hh"
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
    using rtt_meshReaders_test::TestHexFormat;
    
    return SP<TestApp>(new TestHexFormat(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_meshReaders_test
{

using std::string;
using rtt_meshReaders::Hex_Format;

TestHexFormat::TestHexFormat(int argc, char *argv[],
					       std::ostream& os_in)

    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestHexFormat" << endl;
}

string TestHexFormat::version() const
{
    return rtt_meshReaders::release();
}

/*!
 * \brief Tests the CIC-19 Hex mesh format reader.
 *
 */
string TestHexFormat::runTest()
{
    using rtt_meshReaders::Hex_Format;


    if (passed())
    {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}
} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestHexFormat.cc
//---------------------------------------------------------------------------//
