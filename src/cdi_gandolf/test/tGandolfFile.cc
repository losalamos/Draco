//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/tGandolfFile.cc
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:48:56 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "tGandolfFile.hh"
#include "../GandolfFile.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

#include <string>
#include <iostream>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame 
{
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in )
	{
	    using rtt_dsxx::SP;
	    using rtt_cdi_gandolf_test::tGandolfFile;
	    
	    return SP<TestApp>( new tGandolfFile( argc, argv, os_in ) );
	}
} // end namespace rtt_UnitTestFrame

// tGandolfFile Stuff
// ------------------------------------------------------------
namespace rtt_cdi_gandolf_test
{

tGandolfFile::tGandolfFile( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tGandolfFile" << std::endl;
}

std::string tGandolfFile::version() const
{
    return rtt_cdi_gandolf::release();
}

std::string tGandolfFile::runTest()
{
    // Gandolf data filename (IPCRESS format required)
    const std::string op_data_file = "Al_BeCu.ipcress";
    
    // Material identifier.  This data file has two materials: Al and
    // BeCu.  Al has the id tag "10001".
//    const int matid=10001;
    
    // Start the test.

    std::cout << std::endl 
	 << "Testing the GandolfFile component of the cdi_gandolf package." 
	 << std::endl;

    // Create a GandolfFile Object

    std::cout << "Creating a Gandolf File object" << std::endl;
    
    rtt_dsxx::SP<rtt_cdi_gandolf::GandolfFile> spGF;
    spGF = new rtt_cdi_gandolf::GandolfFile( op_data_file );

    // Test the new object to verify the constructor and accessors.

    if ( spGF->getDataFilename() == op_data_file )
	pass() << "Data filename set and retrieved correctly.";
    else
	fail() << "Data filename either not set correctly or not retrieved correctly.";

    if ( spGF->getNumMaterials() == 2 )
	pass() << "Found the correct number of materials in the data file.";
    else
	fail() << "Did not find the correct number of materials in the data file.";
    
    std::cout << std::endl 
	      << "Materials found in the data file:" << std::endl;
    for ( int i=0; i<spGF->getNumMaterials(); ++i )
	std::cout << "  Material " << i << " has the identification number " 
		  << spGF->getMatIDs()[i] << std::endl;

    // Print the test result.
    if (passed()) {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";

} // end of runTest();

} // end namespace rtt_cdi_gandolf_test


//---------------------------------------------------------------------------//
//                         end of tGandolfFile.cc
//---------------------------------------------------------------------------//
