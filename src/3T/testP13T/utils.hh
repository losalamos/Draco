//----------------------------------*-C++-*----------------------------------//
// utils.hh
// Randy M. Roberts
// Wed Nov  4 09:58:01 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_utils_hh__
#define __3T_testP13T_utils_hh__

#include "3T/testP13T/faces.hh"
#include "ds++/SP.hh"

namespace rtt_3T_testP13T
{
 using dsxx::SP;
 
 template<class MT>
 bool isContinuous(const typename MT::fcdsf &rhs,
		   const SP<MT> &spMesh);

 template<class FT>
 double sum(const FT &rhs);
 
 template<class FT>
 typename FT::value_type min(const FT &rhs);
 
 template<class FT>
 typename FT::value_type max(const FT &rhs);

 template<class FT, class OP>
 typename FT::value_type min(const FT &rhs, OP &op);

 template<class FT, class OP>
 typename FT::value_type max(const FT &rhs, OP &op);

} // end namespace rtt_3T_testP13T

// #include "utils.t.cc"

#endif                          // __3T_testP13T_utils_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/utils.hh
//---------------------------------------------------------------------------//
