//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_Traits.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:37:29 2002
 * \brief  Traits for C4 intrinsic types.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_C4_Traits_hh__
#define __c4_C4_Traits_hh__

#include "C4_Tags.hh"

namespace rtt_c4
{
 
//===========================================================================//
/*!
 * \struct C4_Traits
 *
 * This struct and its specializations are used to implement the type-safe
 * default message tags in C4.  Any other type-determined property needed in
 * C4 would also go here.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class T>
struct C4_Traits
{
};

//---------------------------------------------------------------------------//
// SPECIALIZATION OF INTRINSIC ELEMENTAL TYPES
//---------------------------------------------------------------------------//

template<> 
class C4_Traits<char> 
{
    static const int tag = C4_char_Tag;
};

template<> 
class C4_Traits<int> 
{
    static const int tag = C4_int_Tag;
};

template<> 
class C4_Traits<float> 
{
    static const int tag = C4_float_Tag;
};

template<> 
class C4_Traits<double> 
{
    static const int tag = C4_double_Tag;
};

//---------------------------------------------------------------------------//
// SPECIALIZATION OF INTRINSIC POINTER TYPES
//---------------------------------------------------------------------------//

template<> 
class C4_Traits<char *>
{
    static const int tag = C4_char_ptr_Tag;
};

template<>
class C4_Traits<int *> 
{
    static const int tag = C4_int_ptr_Tag;
};

template<> 
class C4_Traits<float *>
{
    static const int tag = C4_float_ptr_Tag;
};

template<> 
class C4_Traits<double *>
{
    static const int tag = C4_double_ptr_Tag;
};

} // end namespace rtt_c4

#endif                          // __c4_C4_Traits_hh__

//---------------------------------------------------------------------------//
//                              end of c4/C4_Traits.hh
//---------------------------------------------------------------------------//
