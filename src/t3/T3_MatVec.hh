//----------------------------------*-C++-*----------------------------------//
// T3_MatVec.hh
// Dave Nystrom
// Thu Oct  2 14:21:13 MDT 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __t3_T3_MatVec_hh__
#define __t3_T3_MatVec_hh__

#include "linalg/PCG_MatVec.hh"

#include "ds++/Mat.hh"

#include "c4/NodeInfo.hh"

class Test_Prob;

//===========================================================================//
// class T3_MatVec - 

// 
//===========================================================================//

template<class T>
class T3_MatVec : public PCG_MatVec<T>, private C4::NodeInfo {

    Test_Prob *prob;

    Mat1<int> ncps, goffs;

    int its;

  public:
    T3_MatVec( Test_Prob *p );
    ~T3_MatVec();

// Get A somehow else.  Then compute b = A x.
    void MatVec( Mat1<T>& b, Mat1<T>& x );

    int get_iterations() const { return its; }
};

#endif                          // __t3_T3_MatVec_hh__

//---------------------------------------------------------------------------//
//                              end of t3/T3_MatVec.hh
//---------------------------------------------------------------------------//
