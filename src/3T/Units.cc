//----------------------------------*-C++-*----------------------------------//
// Units.cc
// Randy M. Roberts
// Tue Mar 17 14:52:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/Units.hh"
#include <limits>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
double Units::minConversion()
{
    static double minConv = std::numeric_limits<double>::min();
    return minConv;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of Units.cc
//---------------------------------------------------------------------------//
