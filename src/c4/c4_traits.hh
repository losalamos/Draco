//----------------------------------*-C++-*----------------------------------//
// c4_traits.hh
// Geoffrey Furnish
// Fri Oct  3 10:01:38 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __c4_c4_traits_hh__
#define __c4_c4_traits_hh__

#include "tags.hh"

//===========================================================================//
// class c4_traits - 

// 
//===========================================================================//

template<class T>
class c4_traits {
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
