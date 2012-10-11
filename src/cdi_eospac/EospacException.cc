//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/EospacException.cc
 * \author Kelly Thompson
 * \date   Fri Apr  6 13:59:06 2001
 * \brief  Implementation file for the cdi_eospac exception handler class.
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "EospacException.hh"
#include <algorithm>

namespace rtt_cdi_eospac
{

    // --------------- //
    // EospacException //
    // --------------- //

    EospacException::EospacException( const std::string& in_what_arg )
	throw() 
	{
	    // allocate space for the error message.
	    int len = in_what_arg.length();
	    message = new char [ len+1 ];
	    
	    // copy the string message to the char* data member.
	    std::copy( in_what_arg.begin(), in_what_arg.end(), 
		       message );

	    // add a terminating character.
	    message[ len ] = 0;

	    // EospacException inherits from std::exception.  This
	    // means that a std::exception is also thrown.
	}

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
    
    // --------------------- //
    // EospacUnknownDataType //
    // --------------------- //

    EospacUnknownDataType::EospacUnknownDataType( 
	const std::string& what_arg ) throw()
	: EospacException( what_arg )
	{
	    // empty
	}

    EospacUnknownDataType::~EospacUnknownDataType()
	throw()
	{
	    // empty
	}

} // end namespace rtt_cdi_eospac

//---------------------------------------------------------------------------//
//                              end of EospacException.cc
//---------------------------------------------------------------------------//
