//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Assert.cc
 * \author Geoffrey Furnish
 * \date   Fri Jul 25 08:41:38 1997
 * \brief  Helper functions for the Assert facility.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <sstream>
#include "Assert.hh"

namespace rtt_dsxx
{

//===========================================================================//
// ASSERTION CLASS MEMBERS
//===========================================================================//

//---------------------------------------------------------------------------//
// Build the error string (PRIVATE)
//---------------------------------------------------------------------------//
std::string assertion::build_message( std::string const &cond, 
				      std::string const &file, 
				      int         const line ) const
{
    using std::string;

    std::ostringstream myMessage;
    myMessage << "Assertion: "
	      << string( cond )
	      << ", failed in "
	      << string( file )
	      << ", line "
	      << line
	      << "." << std::endl;
    return myMessage.str();
}

//===========================================================================//
// FREE FUNCTIONS
//===========================================================================//
/*!
 * \brief Throw a rtt_dsxx::assertion for Require, Check, Ensure macros.
 */
void toss_cookies( std::string const cond, 
		   std::string const file, 
		   int         const line )
{
    throw assertion( cond, file, line );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Throw a rtt_dsxx::assertion for Insist macros.
 */
void insist( std::string const cond, 
	     std::string const msg, 
	     std::string const file, 
	     int         const line)
{
    throw assertion( msg );
}

} // end of rtt_dsxx

//---------------------------------------------------------------------------//
//                              end of Assert.cc
//---------------------------------------------------------------------------//
