//----------------------------------*-C++-*----------------------------------//
// Global.hh
// Thomas M. Evans
// Fri Jan 30 16:23:53 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __imctest_Global_hh__
#define __imctest_Global_hh__

//===========================================================================//
// class Global - 
//
// Date created : 1-10-97
// Purpose      : provide Globally used constants and classes
//
// revision history:
// -----------------
// 0) original
// 1) 2-4-98: added huge and epsilon constants
// 
//===========================================================================//

#include "Names.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <iomanip>

GLOBALSPACE

using std::ostream;
using std::vector;

// CONSTANTS
const double pi = 2.0 * std::asin(1.0);
const double huge = DBL_MAX;
const double epsilon = DBL_EPSILON;

// OVERLOADED OPERATORS
// overloaded operator for printing a vector templated on vector type (VT)
template<class VT>
inline ostream& operator<<(ostream &output, const vector<VT> &object)
{
    using std::setw;
    using std::endl;

    int size = object.size();
    for (int i = 0; i < size; i++)
	output << setw(5) << i << setw(10) << object[i] << endl;
    return output;
}

CSPACE

#endif                          // __imctest_Global_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Global.hh
//---------------------------------------------------------------------------//
