//----------------------------------*-C++-*----------------------------------//
// Global.hh
// Thomas M. Evans
// Mon Jun  8 11:29:48 1998
//---------------------------------------------------------------------------//
// @> Global namespace declarations
//---------------------------------------------------------------------------//

#ifndef __imc_Global_hh__
#define __imc_Global_hh__

//===========================================================================//
// namespace Global - 
//
// Purpose : groups all functions, variables, and constants that are defined 
//           in the IMC::Global namespace into one header.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

#include "mc/Math.hh"
#include "mc/Constants.hh"

namespace rtt_imc 
{
namespace global
{

//---------------------------------------------------------------------------//
// PROBLEM VARIABLES
//---------------------------------------------------------------------------//

//===========================================================================//
/*!
 * \class Type_Switch
 * 
 * \brief Class used to call functions based on template type arguments.
 *
 * This class can be used to simulate partial template specialization for
 * member functions in a class template.  See Alexandrescu, 2001, "Modern C++
 * Design", chapter 2.
 *
 * Note, the Alexandrescu trick works under automatic instantiation rules
 * because functions that are not needed are not instantiated.  However,
 * because we use explicit instantiation, an additional trick is required.
 * The additional trick is to add a "faux" template argument to the function
 * to prevent it from being instantiated as part of the explicit class
 * instantiation.  We will write this up in detail at a later time.
 */
//===========================================================================//

template<class T>
struct Type_Switch
{
    typedef T Type;
};

} // end namespace global
} // end namespace rtt_imc

#endif                          // __imc_Global_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Global.hh
//---------------------------------------------------------------------------//
