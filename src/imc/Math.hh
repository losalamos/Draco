//----------------------------------*-C++-*----------------------------------//
// Math.hh
// Thomas M. Evans
// Fri Mar 20 13:20:50 1998
//---------------------------------------------------------------------------//
// @> IMC::Global namespace math functions
//---------------------------------------------------------------------------//

#ifndef __imctest_Math_hh__
#define __imctest_Math_hh__

//===========================================================================//
// namespace Math - 
//
// Purpose : holds useful mathematical functions
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imctest/Names.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <cassert>

IMCSPACE
GLOBALSPACE

using std::vector;
using std::ostream;

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
	output << setw(5) << i << setw(10) << object[i] << endl;
    return output;
}

//---------------------------------------------------------------------------//
// VECTOR FUNCTIONS
//---------------------------------------------------------------------------//
// do the DOT product between two vectors
template<class VT>
inline VT dot(vector<VT> A, vector<VT> B)
{
    assert (A.size() == B.size());
    VT value = 0.0;
    for (int i = 0; i < A.size(); i++)
	value += A[i] * B[i];
    return value;
}

// do the CROSS product between two vectors
template<class VT>
inline vector<VT> cross(vector<VT> A, vector<VT> B)
{
    assert (A.size() == B.size());
    assert (A.size() == 3);
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
    T return_value = A;
    if (B <= return_value)
	return_value = B;
    return return_value;
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

CSPACE
CSPACE

#endif                          // __imctest_Math_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Math.hh
//---------------------------------------------------------------------------//
