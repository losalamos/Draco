//----------------------------------*-C++-*----------------------------------//
// PreCond.hh
// Geoffrey M. Furnish
// Wed Nov 26 16:37:05 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_test_PreCond_hh__
#define __3T_test_PreCond_hh__

#include "linalg/PCG_PreCond.hh"

#include "ds++/Mat.hh"

//===========================================================================//
// class PreCond - 

// 
//===========================================================================//

template<class Solver>
class PreCond : public PCG_PreCond<typename Solver::NumT>
{
    int ncp;
    int method;
    Solver *solver;

    typedef typename Solver::NumT T;

  public:
    PreCond( int ncp_, int m, Solver *p ) : ncp(ncp_), method(m), solver(p) {}
    void  Left_PreCond( Mat1<T>& x, Mat1<T>&b );
    void Right_PreCond( Mat1<T>& x, Mat1<T>&b );
};

#endif                          // __3T_test_PreCond_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/PreCond.hh
//---------------------------------------------------------------------------//
