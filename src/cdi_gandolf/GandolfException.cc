//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfException.cc
 * \author Kelly Thompson
 * \date   Tue Sep  5 10:47:29 2000
 * \brief  GandolfException class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfException.hh"

#include <sstream> // define std::endl and std::ostringstream

namespace rtt_cdi_gandolf
{

// This function definition allows the user to catch a
// GandolfException exception instead of an exception for a specific
// Gandolf function.  

std::string GandolfException::errorMessage() const
    {
	std::string s1;
	if ( gandolfFunctionName == "gkeys" )
	    {
		gkeysException GandError( errorCode );
		return GandError.errorMessage();
	    }
	else if ( gandolfFunctionName == "gchgrids" )
	    {
		gchgridsException GandError ( errorCode );
		return GandError.errorMessage();
	    }
	else if ( gandolfFunctionName == "ggetgray" )
	    {
		ggetgrayException GandError ( errorCode );
		return GandError.errorMessage();
	    }
	else if ( gandolfFunctionName == "ggetmg" )
	    {
		ggetmgException GandError ( errorCode );
		return GandError.errorMessage();
	    }
	else if ( gandolfFunctionName == "gmatids" )
	    {
		gmatidsException GandError ( errorCode );
		return GandError.errorMessage();
	    }
	else
	    {	
		s1 = "An unknown error has occured in the Gandolf libraries.";
		return s1;
	    }
    } // end of GandolfException::errorMessage();

std::string GandolfException::errorSummary() const
    {
	std::ostringstream outputString;
	outputString << "The Gandolf function named \""
		     << getGandolfFunctionName() << "\""
		     << " returned the error code \""
		     << getErrorCode() << "\"." << std::endl << "\t"
		     << "The message associated with this error code is: " 
		     << std::endl << "\t   "
		     << "\"" << errorMessage() << "\"" 
		     << std::endl;
	return outputString.str();
    }
    
// --------------------- //
// gkeys Exception Class //
// --------------------- //

gkeysException::gkeysException( int errorCode )
    : GandolfException( "gkeys", errorCode)
    {
	// empty
    }

std::string gkeysException::errorMessage() const
    {
	switch ( errorCode ) {
 	case 0: // no errors
 	    return "No error was reported by Gandolf.";
	case -1:
	    return "The requested material ID was not found in the list of material IDs associated with the data file.";
	case -2:
	    return "The requested data key was not found in the list of available keys for this material.";
 	case 1: // IPCRESS file not found.
 	    return "The IPCRESS file was not found.";
	case 2: // File is not IPCRESS.
	    return "The file does not appear to be in IPCRESS format";
	case 3: // Problem reading file
	    return "Having trouble reading the IPCRESS file.";
	case 4: // No keys found for this material.
	    return "No keys were found for this material";
	case 5: // Too many keys found.
	    return "Too many keys for array ( nkeys > kkeys ).";
	default: // unknown error.
	    return "Unknown error returned from Gandolf::gkeys().";
	}
    }

// ------------------------ //
// gchgrids Exception Class //
// ------------------------ //

gchgridsException::gchgridsException( int errorCode )
    : GandolfException( "gchgrids", errorCode)
    {
	// empty
    }

std::string gchgridsException::errorMessage() const
    {
	switch ( errorCode ) {
 	case 0: // no errors
 	    return "No error was reported by Gandolf.";
	case -1: // return with etas, not densities.
	    return "IPCRESS file returned ETAs not densities.";
 	case 1: // IPCRESS file not found.
 	    return "The IPCRESS file was not found.";
	case 2: // File is not IPCRESS.
	    return "The file does not appear to be in IPCRESS format";
	case 3: // Problem reading file
	    return "Having trouble reading the IPCRESS file.";
	case 4: // Inconsistent gray grids, mg not checked
	    return "Gray grid inconsistent with the temp/density grid.";
	case 5: // ngray != nt*nrho, mg not checked
	    return "Wrong number of gray opacities found (ngray != nt*nrho)." ;
	case 6: // inconsistent mg grid.
	    return "MG grid inconsistent with the temp/density/hnu grid.";
	case 7: //  nmg != nt*nrho*(nhnu-1).
	    return "Wrong number of MG opacities found (nmg != nt*nrho*(nhnu-1)).";
	default: // unknown error.
	    return "Unknown error returned from Gandolf::gchgrids().";
	}
    }


// ------------------------ //
// ggetgray Exception Class //
// ------------------------ //

ggetgrayException::ggetgrayException( int errorCode )
    : GandolfException( "ggetgray", errorCode)
    {
	// empty
    }

std::string ggetgrayException::errorMessage() const
    {
	switch ( errorCode ) {
 	case 0: // no errors
 	    return "No error was reported by Gandolf.";
	case -1: // return with etas, not densities.
	    return "IPCRESS file returned ETAs not densities.";
 	case 1: // IPCRESS file not found.
 	    return "The IPCRESS file was not found.";
	case 2: // File is not IPCRESS.
	    return "The file does not appear to be in IPCRESS format";
	case 3: // Problem reading file
	    return "Having trouble reading the IPCRESS file.";
	case 4: // Data not found
	    return "Requested data not found.  Check nt, nrho, ngray.";
	case 5: // Data larger than allocated arrays.
	    return "Data found is larger than allocated array size.";
	case 6: // Data size not equal to nt*nrho
	    return "Data size not equal to expected size (ndata != nt*nrho)";
	case 7: // Opacity requested but no table loaded.
	    return "The gray opacity data table is not currently available.";
	default: // unknown error.
	    return "Unknown error returned from Gandolf::ggetgray().";
	}
    }

// ---------------------- //
// ggetmg Exception Class //
// ---------------------- //

ggetmgException::ggetmgException( int errorCode )
    : GandolfException( "ggetmg", errorCode)
    {
	// empty
    }

std::string ggetmgException::errorMessage() const
    {
	switch ( errorCode ) {
 	case 0: // no errors
 	    return "No error was reported by Gandolf.";
	case -1: // return with etas, not densities.
	    return "IPCRESS file returned ETAs not densities.";
 	case 1: // IPCRESS file not found.
 	    return "The IPCRESS file was not found.";
	case 2: // File is not IPCRESS.
	    return "The file does not appear to be in IPCRESS format";
	case 3: // Problem reading file
	    return "Having trouble reading the IPCRESS file.";
	case 4: // Data not found
	    return "Requested data not found.  Check nt, nrho, ngray.";
	case 5: // Data larger than allocated arrays.
	    return "Data found is larger than allocated array size.";
	case 6: // Data size not equal to nt*nrho
	    return "Data size not equal to expected size (ndata != nt*nrho*(nhnu-1))";
	default: // unknown error.
	    return "Unknown error returned from Gandolf::ggetmg().";
	}
    }

// ----------------------- //
// gmatids Exception Class //
// ----------------------- //

gmatidsException::gmatidsException( int errorCode )
    : GandolfException( "gmatids", errorCode)
    {
	// empty
    }

std::string gmatidsException::errorMessage() const
    {
	switch ( errorCode ) {
 	case 0: // no errors
 	    return "No error was reported by Gandolf.";
	case -1: // Filename is too long.
	    return "The filename given to Gandolf has too many characters (maxlen=80).";
 	case 1: // IPCRESS file not found.
 	    return "The IPCRESS file was not found.";
	case 2: // File is not IPCRESS.
	    return "The file does not appear to be in IPCRESS format";
	case 3: // Problem reading file
	    return "Having trouble reading the IPCRESS file.";
	case 4: // No material ID's found in file.
	    return "No material ID's were found in the IPCRESS data file.";
	case 5: // too many matids found ( nmat > kmat )
	    return "Too many materials were found in the data file ( nmat > kmat ).";
	default: // unknown error.
	    return "Unknown error returned from Gandolf::gmatids().";
	}
    }

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
// end of GandolfException.cc
//---------------------------------------------------------------------------//
