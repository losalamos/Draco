//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/tGandolfFile.cc
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:39:39 2001
 * \brief  Gandolf file test
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_gandolf_test.hh"
#include "ds++/Release.hh"
#include "../GandolfFile.hh"
#include "../GandolfException.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_cdi_gandolf::GandolfFile;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
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
 *
 * 2010-12-01: When using Draco Modules (draco/environment/Modules), you
 * simply need to load the gandolf module and the build system will
 * automatically find libgandolf by looking at the environment variable
 * GANDOLF_LIB_DIR. 
 */
void gandolf_file_test()
{
    // Gandolf data filename (IPCRESS format required)
    const std::string op_data_file = "Al_BeCu.ipcress";
	    
    // Start the test.
	    
    std::cout << std::endl 
	      << "Testing the GandolfFile component of the "
	      << "cdi_gandolf package." << std::endl;
	    
    // Create a GandolfFile Object
	    
    std::cout << "Creating a Gandolf File object" << std::endl;
	    
    SP<GandolfFile> spGF;
    try
    {
	spGF = new rtt_cdi_gandolf::GandolfFile( op_data_file );
    }
    catch ( const rtt_cdi_gandolf::gmatidsException& GandError )
    {
	FAILMSG(GandError.what());
	FAILMSG(GandError.errorSummary());
	return;
    }

    // Test the new object to verify the constructor and accessors.
	    
    std::vector<int> matIDs = spGF->getMatIDs();
    if ( matIDs[0] == 10001 && matIDs[1] == 10002 )
    {
	PASSMSG("Found two materials in IPCRESS file with expected IDs.");
    }
    else
    {
	FAILMSG("Did not find materials with expected IDs in IPCRESS file.");
    }
	    
    if ( spGF->materialFound( 10001 ) )
    {
	PASSMSG("Looks like material 10001 is in the data file.");
    }
    else
    {
	FAILMSG("Can't find material 10001 in the data file.");
    }
	    
    if ( spGF->materialFound( 5500 ) ) // should fail
    {
	ostringstream message;
	message << "Material 5500 shouldn't exist in the data file." 
		<< "\n\tLooks like we have a problem.";
	FAILMSG(message.str());
    }
    else
    {
	ostringstream message; 
	message << "Access function correctly identified material 5500"
		<< "\n\tas being absent from IPCRESS file.";
	PASSMSG(message.str());
    }
	    
    if ( spGF->getDataFilename() == op_data_file )
    {
	PASSMSG("Data filename set and retrieved correctly.");
    }
    else
    {
	ostringstream message;
	message << "Data filename either not set correctly or not "
		<< "retrieved correctly.";
	FAILMSG(message.str());
    }
	    
    if ( spGF->getNumMaterials() == 2 )
    {
	PASSMSG("Found the correct number of materials in the data file.");
    }
    else
    {
	ostringstream message;
	message << "Did not find the correct number of materials in "
		<< "the data file.";
	FAILMSG(message.str());
    }
	    
    std::cout << std::endl 
	      << "Materials found in the data file:" << std::endl;
	    
    for ( size_t i=0; i<spGF->getNumMaterials(); ++i )
	std::cout << "  Material " << i << " has the identification number " 
		  << spGF->getMatIDs()[i] << std::endl;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	gandolf_file_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tGandolfFile, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout <<     "\n*********************************************";
    if (rtt_cdi_gandolf_test::passed) 
        cout << "\n**** tGandolfFile Test: PASSED"; 
    cout <<     "\n*********************************************\n";
    cout <<     "\nDone testing tGandolfFile." << endl;

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tGandolfFile.cc
//---------------------------------------------------------------------------//
