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
std::string assertion::build_message( std::string const & cond, 
				      std::string const & file, 
				      int         const line ) const
{
    std::ostringstream myMessage;
    myMessage << "Assertion: "
	      << cond
	      << ", failed in "
	      << file
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
void toss_cookies( std::string const & cond, 
		   std::string const & file, 
		   int         const line )
{
    throw assertion( cond, file, line );
}


void 
toss_cookies_ptr(char const * const cond, 
		 char const * const file, 
		 int  const line )
{
    throw assertion( cond, file, line );
}


//---------------------------------------------------------------------------//
/*! 
 * \brief Throw a rtt_dsxx::assertion for Insist macros.
 */
void insist( std::string const & cond, 
	     std::string const & msg, 
	     std::string const & file, 
	     int         const line)
{
    std::ostringstream myMessage;
    myMessage <<  "Insist: " << cond << ", failed in "
	      << file << ", line " << line << "." << std::endl
	      << "The following message was provided:" << std::endl
	      << "\"" << msg << "\"" << std::endl;
    throw assertion( myMessage.str() );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Throw a rtt_dsxx::assertion for Insist_ptr macros.
 * Having a (non-inlined) version that takes pointers prevents the compiler
 * from having to construct std::strings from the pointers each time.  This
 * is particularly important for things like rtt_dsxx::SP::operator->, that
 * (a) have an insist in them, (b) don't need complicated strings and (c) are
 * called frequently.
 */
void insist_ptr( char const * const cond, 
		 char const * const msg, 
		 char const * const file, 
		 int          const line)
{
    // Call the other insist for consistency
    insist(cond, msg, file, line);
}


} // end of rtt_dsxx

//---------------------------------------------------------------------------//
//                              end of Assert.cc
//---------------------------------------------------------------------------//
