//----------------------------------*-C++-*----------------------------------//
// Source_Init.hh
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Source_Init_hh__
#define __imctest_Source_Init_hh__

//===========================================================================//
// class Source_Init - 
//
// Purpose : do source initialization for the source builder and the parallel 
//           plan
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imctest/Names.hh"

template<class MT>
class Source_Init
{
private:

public:
    Source_Init(SP<MT> mesh);

  // source initialyzer function
    void Initialize();

#endif                          // __imctest_Source_Init_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Source_Init.hh
//---------------------------------------------------------------------------//
