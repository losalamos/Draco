//----------------------------------*-C++-*----------------------------------//
// DbgOstream.hh
// Randy M. Roberts
// Wed Jun 23 14:35:32 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_DbgOstream_hh__
#define __matprops_DbgOstream_hh__

#include <iostream>
#include "NullOstream.hh"

namespace rtt_matprops
{
 
//===========================================================================//
// class DbgOstream - 
//
// Purpose : Provide a sometimes on, sometimes off ostream
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class DbgOstream 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef std::basic_ostream<char, std::char_traits<char> > ostream_type;

    // DATA
    
    std::ostream *os;
    bool isNull;

  public:

    // CREATORS
    
    DbgOstream(std::ostream &os_)
	: os(&os_), isNull(false)
    {
	// empty
    }
    
    DbgOstream()
	: os(new NullOstream), isNull(true)
    {
	// empty
    }

    ~DbgOstream()
    {
	if (isNull)
	    delete os;
    }

    // MANIPULATORS
    
    template<class T>
    DbgOstream &operator<<(const T& t)
    {
	*os << t;
	return *this;
    }

    DbgOstream &operator<<(ostream_type& (*pf) (ostream_type&))
    {
	*os << pf;
	return *this;
    }

    // ACCESSORS

  private:

    // DISALLOWED OPERATIONS
    
    DbgOstream(const DbgOstream &rhs);
    DbgOstream& operator=(const DbgOstream &rhs);

    // IMPLEMENTATION
};

} // end namespace rtt_matprops

#endif                          // __matprops_DbgOstream_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/DbgOstream.hh
//---------------------------------------------------------------------------//
