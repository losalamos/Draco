//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/tGandolfFile.cc
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:48:56 2000
 * \brief  Implementation file for tGandolfFile
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

// cdi_gandolf dependencies
#include "tGandolfFile.hh"
#include "../GandolfFile.hh"
#include "../GandolfException.hh"
#include "../Release.hh"

// Draco dependencies
#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

// STL dependencies
#include <string>
#include <iostream>

// DEBUG dependencies
// #include <iomanip>

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
    
    //===========================================================================
    /*!
     * \brief Tests the GandolfFile constructor and access routines.
     *
     * Information about the GANDOLF IPCRESS Interface Library may be
     * found on X-Division's Physical Data Team web site:
     *
     * http://velvet.lanl.gov/PROJECTS/DATA/atomic/gandolf/intro.lasso
     *
     * We have a slightly modified copy of the gandolf libraries
     * (explicitly added compiler libraries) located at:
     *
     * /radtran/vendors/gandolf
     *
     * To use this package Draco must be compiled with the one of the
     * following tags (depending on the machine architecture):
     *
     * --with-gandolf-lib=/radtran/vendors/gandolf/IRIX64/lib64
     * --with-gandolf-lib=/radtran/vendors/gandolf/IRIX64/lib32
     * --with-gandolf-lib=/radtran/vendors/gandolf/Linux
     * --with-gandolf-lib=/radtran/vendors/gandolf/SunOS
     */
    //===========================================================================
    std::string tGandolfFile::runTest()
	{
	    // Gandolf data filename (IPCRESS format required)
	    const std::string op_data_file = "Al_BeCu.ipcress";
	    
	    // Start the test.
	    
	    std::cout << std::endl 
		      << "Testing the GandolfFile component of the "
		      << "cdi_gandolf package." << std::endl;
	    
	    // Create a GandolfFile Object
	    
	    std::cout << "Creating a Gandolf File object" << std::endl;
	    
	    rtt_dsxx::SP<rtt_cdi_gandolf::GandolfFile> spGF;
	    try
		{
		    spGF = new rtt_cdi_gandolf::GandolfFile( op_data_file );
		}
	    catch ( const rtt_cdi_gandolf::gmatidsException& GandError )
		{
		    fail() << std::endl << "\t" << GandError.what();
		    fail() << std::endl << "\t" << GandError.errorSummary();
		    return "Some tests failed.";
		}
	    
	    // Test the new object to verify the constructor and accessors.
	    
	    std::vector<int> matIDs = spGF->getMatIDs();
	    if ( matIDs[0] == 10001 && matIDs[1] == 10002 )
		pass () << "Found two materials in IPCRESS file with expected IDs.";
	    else
		fail () << "Did not find materials with expected IDs in IPCRESS file.";
	    
	    if ( spGF->materialFound( 10001 ) )
		pass () << "Looks like material 10001 is in the data file.";
	    else
		fail() << "Can't find material 10001 in the data file.";
	    
	    if ( spGF->materialFound( 5500 ) ) // should fail
		fail() << "Material 5500 shouldn't exist in the data file." 
		       << "\n\tLooks like we have a problem.";
	    else
		pass() << "Access function correctly identified material 5500"
		       << "\n\tas being absent from IPCRESS file.";
	    
	    if ( spGF->getDataFilename() == op_data_file )
		pass() << "Data filename set and retrieved correctly.";
	    else
		fail() << "Data filename either not set correctly or not "
		       << "retrieved correctly.";
	    
	    if ( spGF->getNumMaterials() == 2 )
		pass() << "Found the correct number of materials in the data file.";
	    else
		fail() << "Did not find the correct number of materials in "
		       << "the data file.";
	    
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
