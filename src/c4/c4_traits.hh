//----------------------------------*-C++-*----------------------------------//
// c4_traits.hh
// Geoffrey Furnish
// Fri Oct  3 10:01:38 1997
//---------------------------------------------------------------------------//
// @> Traits specializations for intrinsic types in C4.
//---------------------------------------------------------------------------//

#ifndef __c4_c4_traits_hh__
#define __c4_c4_traits_hh__

#include "tags.hh"

//===========================================================================//
// class c4_traits - Define properties of types for C4

// This class and its specializations are used to implement the type-safe
// default message tags in C4.  Any other type-determined property needed in
// C4 would also go here.
//===========================================================================//

template<class T>
class c4_traits {
};

// Intrinsic elemental types.

template<> class c4_traits<char> {
  public:
    static const int Tag = C4_char_Tag;
};

template<> class c4_traits<int> {
  public:
    static const int Tag = C4_int_Tag;
};

template<> class c4_traits<float> {
  public:
    static const int Tag = C4_float_Tag;
};

template<> class c4_traits<double> {
  public:
    static const int Tag = C4_double_Tag;
};

// Intrinsic pointer types.

template<> class c4_traits<char *> {
  public:
    static const int Tag = C4_char_ptr_Tag;
};

template<> class c4_traits<int *> {
  public:
    static const int Tag = C4_int_ptr_Tag;
};

template<> class c4_traits<float *> {
  public:
    static const int Tag = C4_float_ptr_Tag;
};

template<> class c4_traits<double *> {
  public:
    static const int Tag = C4_double_ptr_Tag;
};

#endif                          // __c4_c4_traits_hh__

//---------------------------------------------------------------------------//
//                              end of c4/c4_traits.hh
//---------------------------------------------------------------------------//
