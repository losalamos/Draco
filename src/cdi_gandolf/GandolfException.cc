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
    
    //========================================================================
    /*!
     * \class GandolfException
     *
     * \brief This class handles exceptions thrown by GandolfOpacity when
     *        it calls the Gandolf library functions.
     *
     * \sa GandolfException.hh for more details
     */
    //========================================================================

    GandolfException::GandolfException( 
	const std::string& in_what_arg, 
	std::string in_gandolfFunctionName,
	int in_errorCode ) throw()

	: gandolfFunctionName ( in_gandolfFunctionName ),
	errorCode ( in_errorCode ),
	std::exception ()
// This works for KCC 3.3 on IRIX64
// 	std::exception ( "Gandolf Exception: " + in_what_arg )

	{
	    // empty
	}

    GandolfException::~GandolfException( ) throw()
	{
	    // message is allocated by what().  If what() was never
	    // called then message is a null pointer.
	    if ( message )
		{
		    delete [] message;
		    message = NULL;
		}
	}

    const char* GandolfException::what() const throw() 
	{
	    // if message has not been initalized then we need to call 
	    // errorSummary() to construct the text string.  This
	    // string is then converted into a char* and cached as a
	    // data member of GandolfException.  Memory deallocation
	    // is handled by ~GandolfException().
	    
	    if ( ! message )
		{
		    // access the string.
		    std::string msg = errorSummary();
		    int len = msg.length();
		    // allocate space for converting the string to a char*
		    message = new char [ len + 1 ];
		    msg.copy( message, len, 0 );
		    // add a string termination character.
		    message[ len ] = 0;
		}

	    return message; 
	}

    std::string GandolfException::errorSummary() const
	{
	    std::ostringstream outputString;
	    outputString << "The Gandolf function named \""
			 << getGandolfFunctionName() << "\""
			 << " returned the error code \""
			 << getErrorCode() << "\"." << std::endl << "\t"
			 << "The message associated with this error code is: " 
			 << std::endl << "\t   "
			 << "\"" << getErrorMessage() << "\"" 
			 << std::endl;
 	    return outputString.str();
	}

    std::string GandolfException::getErrorMessage() const
	{
	    // This should never be accessed.  The derived exception
	    // objects should override this virtual function.
 	    return "bogus error message.";
	}


    //========================================================================
    /*!
     * \class gkeysException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function gkeys().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     */
    //========================================================================

    gkeysException::gkeysException( const std::string& what_arg )
	: GandolfException ( what_arg, "gkeys()" )
	{
	    // empty
	}

    gkeysException::gkeysException( int gkeysError )
	: GandolfException ( "gkeysError", "gkeys()", gkeysError )
	{	   
	    // empty
	}

    gkeysException::~gkeysException() throw()
	{
	    // empty
	}

    std::string gkeysException::getErrorMessage() const
	{
	    // look up errorCode to find associated message.
	    switch ( errorCode ) 
		{
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


   //========================================================================
    /*!
     * \class gchgridsException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function gchgrids().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     */
    //========================================================================

    gchgridsException::gchgridsException( const std::string& what_arg )
	: GandolfException ( what_arg, "gchgrids()" )
	{
	    // empty
	}

    gchgridsException::gchgridsException( int gchgridsError )
	: GandolfException ( "gchgridsError", "gchgrids()", gchgridsError ) 
	{
	    // empty
	}

    gchgridsException::~gchgridsException() throw()
	{
	    // empty
	}

    std::string gchgridsException::getErrorMessage() const
	{
	    // look up errorCode to find associated message.
	    switch ( errorCode ) 
		{
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


    //========================================================================
    /*!
     * \class ggetgrayException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function ggetgray().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     */
    //========================================================================

    ggetgrayException::ggetgrayException( const std::string& what_arg )
	: GandolfException ( what_arg, "ggetgray()" )
	{
	    // empty
	}

    ggetgrayException::ggetgrayException( int ggetgrayError )
	: GandolfException ( "ggetgrayError", "ggetgray()", ggetgrayError )
	{
	    // empty
	}

    ggetgrayException::~ggetgrayException() throw()
	{
	    // empty
	}

    std::string ggetgrayException::getErrorMessage() const
	{
	    // look up errorCode to find associated message.
	    switch ( errorCode ) 
		{
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
    

    //========================================================================
    /*!
     * \class ggetmgException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function ggetmg().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     */
    //========================================================================

    ggetmgException::ggetmgException( const std::string& what_arg )
	: GandolfException ( what_arg, "ggetmg()" )
	{
	    // empty
	}

    ggetmgException::ggetmgException( int ggetmgError )
	: GandolfException( "ggetmgError", "ggetmg()", ggetmgError )
	{
	    // empty
	}

    ggetmgException::~ggetmgException() throw()
	{
	    // empty
	}

    std::string ggetmgException::getErrorMessage() const
	{
	    // look up errorCode to find associated message.
	    switch ( errorCode ) 
		{
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


    //========================================================================
    /*!
     * \class gmatidsException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function gmatids().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     */
    //========================================================================

    gmatidsException::gmatidsException( const std::string& what_arg )
	: GandolfException ( what_arg, "gmatids()" )
	{
	    // empty
	}

    gmatidsException::gmatidsException( int gmatidsError )
 	: GandolfException ( "gmatidsError", "gmatids()", gmatidsError )
	{
	    // empty
	}

    gmatidsException::~gmatidsException() throw()
	{
	    // empty
	}

    std::string gmatidsException::getErrorMessage() const
	{
	    // look up errorCode to find associated message.
	    switch ( errorCode )
		{
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
