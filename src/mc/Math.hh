//----------------------------------*-C++-*----------------------------------//
// Math.hh
// Thomas M. Evans and Todd J. Urbatsch
// Fri Mar 20 13:20:50 1998
//---------------------------------------------------------------------------//
// @> MC::global namespace math functions
//---------------------------------------------------------------------------//

#ifndef __mc_Math_hh__
#define __mc_Math_hh__

//===========================================================================//
// namespace Math - 
//
// Purpose : holds useful mathematical functions
//
// revision history:
// -----------------
// 0) original
// 1) 04-13-99 : moved into mc package 
// 2) 08-04-00 : modified soft_equiv to handle comparisons to zero.
// 3) 08-24-00 : added mod_with_2e9 function.
// 4) 07-31-01 : modified mod_with_2e9 function to integer_modulo with an 
//               argument for the integer to mod against
//
//===========================================================================//

#include "ds++/Assert.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <iomanip>

namespace rtt_mc 
{
namespace global
{

using std::vector;
using std::ostream;

//---------------------------------------------------------------------------//
// CHECK TWO FLOATING POINT NUMBERS FOR EQUALITY
//---------------------------------------------------------------------------//

template<class FPT>
inline bool soft_equiv(const FPT &value, const FPT &reference, 
		       const FPT precision = 1.0e-12)
{
    using std::fabs;
    bool passed = false;

    if (fabs(value - reference) < precision * fabs(reference))
	passed = true;
    else
	passed = false;

    // second chance for passing if reference is identically zero
    if (reference == 0 && !passed)
	if (fabs(value) < precision) passed = true;

    // check for both zeroes
    if (reference == 0 && value == 0 && !passed) passed = true;

    return passed;
}

template<>
inline bool soft_equiv(const int &value, const int &reference,
		       const int precision)
{
    Insist (0, "Can't do a soft compare with integers!");
    return false;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//
// overloaded operator for printing a vector templated on vector type (VT)

template<class VT>
inline ostream& operator<<(ostream &output, const vector<VT> &object)
{
    using std::setw;
    using std::endl;

    int size = object.size();
    for (int i = 0; i < size; i++)
	output << setw(5) << i << "\t" << object[i] << endl;
    return output;
}

//---------------------------------------------------------------------------//
// VECTOR FUNCTIONS
//---------------------------------------------------------------------------//
// do the DOT product between two vectors

template<class VT>
inline VT dot(vector<VT> A, vector<VT> B)
{
    Check (A.size() == B.size());
    VT value = 0.0;
    for (int i = 0; i < A.size(); i++)
	value += A[i] * B[i];
    return value;
}

//---------------------------------------------------------------------------//
// do the CROSS product between two vectors

template<class VT>
inline vector<VT> cross(vector<VT> A, vector<VT> B)
{
    Check (A.size() == B.size());
    Check (A.size() == 3);
    vector<VT> C(3);
    C[0] = A[1] * B[2] - A[2] * B[1];
    C[1] = A[2] * B[0] - A[0] * B[2];
    C[2] = A[0] * B[1] - A[1] * B[0];
    return C;
}

//---------------------------------------------------------------------------//
// GENERAL MATH FUNCTIONS
//---------------------------------------------------------------------------//
// return the minimum value

template<class T>
inline T min(T A, T B)
{
    return A < B ? A : B;
}

template<class T>
inline T min(T A, T B, T C)
{
    T return_value = A;
    if (B <= return_value)
	return_value = B;
    if (C <= return_value)
	return_value = C;
    return return_value;
}

template<class T>
inline T min(T A, T B, T C, T D)
{
    T return_value = A;
    if (B <= return_value)
	return_value = B;
    if (C <= return_value)
	return_value = C;
    if (D <= return_value)
	return_value = D;
    return return_value;
}
	
template<class T>
inline T min(T A, T B, T C, T D, T E)
{
    T return_value = A;
    if (B <= return_value)
	return_value = B;
    if (C <= return_value)
	return_value = C;
    if (D <= return_value)
	return_value = D;
    if (E <= return_value)
	return_value = E;
    return return_value;
}

//---------------------------------------------------------------------------//
// return the max value

template<class T>
inline T max(T A, T B)
{
    return A > B ? A : B;
}

//---------------------------------------------------------------------------//
// modulo an integer with a big integer

inline int integer_modulo(const int value, unsigned int big_int)
{
    Require (value >= 0);

    int mod_value = value % big_int;
    
    Ensure (mod_value < big_int);
    Ensure (mod_value >= 0);

    return mod_value;
}

} // end namespace global
} // end namespace rtt_mc


// Macros to integer_modulo
#define INTEGER_MODULO_1E9(v) rtt_mc::global::integer_modulo(v, 1000000000)
#define INTEGER_MODULO_2E9(v) rtt_mc::global::integer_modulo(v, 2000000000)

#endif                          // __mc_Math_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Math.hh
//---------------------------------------------------------------------------//
