//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfFile.cc
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:15:49 2000
 * \brief  Implementation file for GandolfFile class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfFile.hh"
#include "GandolfWrapper.hh"

namespace rtt_cdi_gandolf
{

    GandolfFile::GandolfFile( const std::string& gandolfDataFilename_ )
	: dataFilename( gandolfDataFilename_ )
	{
	    // Gandolf will only look at the first "maxDataFilenameLength"
	    // characters of the data filename.  We must require that the
	    // given name is less than "maxDataFilenameLength" characters.
	    Require( dataFilename.length() < wrapper::maxDataFilenameLength );
	 
	    // local variables
	    int errorCode = 0;
	    
	    // This call to Gandolf validates the datafile and if
	    // successful returns a list of material identifiers for which 
	    // the materials that exist in the data file.
	    wrapper::gmatids( dataFilename, matIDs, wrapper::maxMaterials,
			      numMaterials, errorCode ); 
	    
	    // Abort if Gandolf returned an error.
	    switch ( errorCode ) {
	    case 0: // no errors
		break;
	    case 1: // IPCRESS file not found.
		Insist( false, "The IPCRESS file was not found.");
		break;
	    case 2: // File is not IPCRESS.
		Insist( false, "The file does not appear to be in IPCRESS format");
		break;
	    case 3: // Problem reading file
		Insist( false, "Having trouble reading the IPCRESS file.");
		break;
	    case 4: // No material ID's found in file.
		Insist( false, "No material ID's were found in the IPCRESS data file.");
		break;
	    case 5: // too many matids found ( nmat > kmat )
		Insist( false, "Too many materials were found in the data file ( nmat > kmat ).");
		break;
	    default: // unknown error.
		Insist( false, "Unknown error returned from Gandolf::gmatids().");
		break;
	    }
   
	}

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
//                              end of GandolfFile.cc
//---------------------------------------------------------------------------//
