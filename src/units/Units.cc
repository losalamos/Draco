//----------------------------------*-C++-*----------------------------------//
// Units.cc
// Randy M. Roberts
// Tue Mar 17 14:52:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Units.hh"
#include <limits>

namespace XTM
{

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

bool operator==(const Units &op1, const Units &op2)
{
    return op1.lengthConversion == op2.lengthConversion &&
		 op1.massConversion == op2.massConversion &&
		 op1.timeConversion == op2.timeConversion &&
		 op1.temperatureConversion == op2.temperatureConversion;
}

} // end namespace XTM

//---------------------------------------------------------------------------//
//                              end of Units.cc
//---------------------------------------------------------------------------//
