//----------------------------------*-C++-*----------------------------------//
// MatVec_3T.hh
// Geoffrey M. Furnish
// Wed Nov 26 16:18:46 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_test_MatVec_3T_hh__
#define __3T_test_MatVec_3T_hh__

#include "ds++/Mat.hh"
#include "ds++/SP.hh"

#include "linalg/PCG_MatVec.hh"

#include "c4/NodeInfo.hh"

class Test_Prob;

//===========================================================================//
// class MatVec_3T - 

// 
//===========================================================================//

template<class MT, class Solver>
class MatVec_3T : public PCG_MatVec<typename Solver::NumT>,
		  private C4::NodeInfo
{
    SP<MT> spm;
    Solver *solver;

    Mat1<int> ncps, goffs;

    int its;

    typedef typename Solver::NumT T;

  public:
    MatVec_3T( const SP<MT>& spm_, Solver *p );
//    ~MatVec_3T();

// Get A somehow else.  Then compute b = A x.
    void MatVec( Mat1<T>& b, const Mat1<T>& x );

    int get_iterations() const { return its; }
};

#endif                          // __3T_test_MatVec_3T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/MatVec_3T.hh
//---------------------------------------------------------------------------//
