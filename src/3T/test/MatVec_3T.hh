//----------------------------------*-C++-*----------------------------------//
// MatVec_3T.hh
// Geoffrey M. Furnish
// Wed Nov 26 16:18:46 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_test_MatVec_3T_hh__
#define __3T_test_MatVec_3T_hh__

#include "linalg/PCG_MatVec.hh"

#include "ds++/Mat.hh"

#include "c4/NodeInfo.hh"

class Test_Prob;

//===========================================================================//
// class MatVec_3T - 

// 
//===========================================================================//

template<class Problem>
class MatVec_3T : public PCG_MatVec<typename Problem::NumT>, private C4::NodeInfo {

    Problem *prob;

    Mat1<int> ncps, goffs;

    int its;

    typedef typename Problem::NumT T;

  public:
    MatVec_3T( Problem *p );
    ~MatVec_3T();

// Get A somehow else.  Then compute b = A x.
    void MatVec( Mat1<T>& b, Mat1<T>& x );

    int get_iterations() const { return its; }
};

#endif                          // __3T_test_MatVec_3T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/MatVec_3T.hh
//---------------------------------------------------------------------------//
