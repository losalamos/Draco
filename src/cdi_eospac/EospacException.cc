//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/EospacException.cc
 * \author Kelly Thompson
 * \date   Fri Apr  6 13:59:06 2001
 * \brief  Implementation file for the cdi_eospac exception handler class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "EospacException.hh"

#include <iostream>

namespace rtt_cdi_eospac
{
    EospacException::EospacException( const std::string& in_what_arg )
	throw() 
	{
	    int len = in_what_arg.length();
	    message = new char [ len+1 ];
	    std::copy( in_what_arg.begin(), in_what_arg.end(), 
		       message );
	    // add a terminating character.
	    message[ len ] = 0;
	}

//     EospacException::EospacException( int errorCode )
// 	throw() 
// 	{
// 	    // All we have is a non-zero error code returned from the
// 	    // EOSPAC library.  We now query EOSPAC for the associated 
// 	    // error message


// 	    int len = in_what_arg.length();
// 	    message = new char [ len+1 ];
// 	    std::copy( in_what_arg.begin(), in_what_arg.end(), 
// 		       message );
// 	    // add a terminating character.
// 	    message[ len ] = 0;
// 	}

    EospacException::~EospacException() throw()
	{
	    // message is allocated by what().  If what() was never
	    // called then message is null.
	    if ( message ) delete [] message;
	}

    const char* EospacException::what() const throw()
	{
	    return message;
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
//                              end of EospacException.cc
//---------------------------------------------------------------------------//
