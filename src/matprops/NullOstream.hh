//----------------------------------*-C++-*----------------------------------//
// NullOstream.hh
// Randy M. Roberts
// Wed Jun 23 14:34:23 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_NullOstream_hh__
#define __matprops_NullOstream_hh__

#include <iostream>

namespace rtt_matprops
{
 
//===========================================================================//
// class NullOstream - 
//
// Purpose : Provide No-Op ostream
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class NullOstream : public std::ostream
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    NullOstream() : std::ostream(0) { };

    // MANIPULATORS
    
    // ACCESSORS

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_matprops

#endif                          // __matprops_NullOstream_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/NullOstream.hh
//---------------------------------------------------------------------------//
