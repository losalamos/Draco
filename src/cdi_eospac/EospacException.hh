//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/EospacException.hh
 * \author Kelly Thompson
 * \date   Fri Apr  6 13:59:06 2001
 * \brief  Header file for the cdi_eospac exception handler class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_EospacException_hh__
#define __cdi_eospac_EospacException_hh__

#include <stdexcept>
#include <string>

namespace rtt_cdi_eospac
{
 
//===========================================================================//
/*!
 * \class EospacException
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

    class EospacException : public std::exception
    {
	
	// NESTED CLASSES AND TYPEDEFS
	
	// DATA

	mutable char* message;
	
      public:
	
	// CREATORS
	
	EospacException( const std::string& what_arg ) throw();
// 	EospacException( int errorCode ) throw();
	// EospacException(const EospacException &rhs);
	virtual ~EospacException() throw();
	
	// MANIPULATORS
	
	// EospacException& operator=(const EospacException &rhs);
	
	// ACCESSORS
	
	virtual const char* what() const throw();

      private:
	
	// IMPLEMENTATION
    };
    
} // end namespace rtt_cdi_eospac

#endif // __cdi_eospac_EospacException_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/EospacException.hh
//---------------------------------------------------------------------------//
