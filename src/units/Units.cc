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
    
//------------------------------------------------------------------------//
// minConversion:
//    A class utility that returns the minimum allowed conversion factor.
//------------------------------------------------------------------------//

double Units::minConversion()
{
    static double minConv = std::numeric_limits<double>::min();
    return minConv;
}

Units operator/(const Units &op1, const Units &op2)
{
    return Units(op1.lengthConversion      / op2.lengthConversion,
		 op1.massConversion        / op2.massConversion,
		 op1.timeConversion        / op2.timeConversion,
		 op1.temperatureConversion / op2.temperatureConversion);
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of Units.cc
//---------------------------------------------------------------------------//
