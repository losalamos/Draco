//----------------------------------*-C++-*----------------------------------//
// MatVec_3T.hh
// Geoffrey M. Furnish
// Wed Nov 26 16:18:46 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_MatVec_3T_hh__
#define __3T_MatVec_3T_hh__

#include "ds++/Mat.hh"
#include "ds++/SP.hh"

#include "linalg/PCG_MatVec.hh"

#include "c4/NodeInfo.hh"

//===========================================================================//
// class MatVec_3T - 

// 
//===========================================================================//

template<class Solver>
class MatVec_3T : public PCG_MatVec<typename Solver::NumT>
{
    Solver *solver;

    int its;

    typedef typename Solver::NumT T;

  public:
    MatVec_3T( Solver *p );

// Get A somehow else.  Then compute b = A x.
    void MatVec( dsxx::Mat1<T>& b, const dsxx::Mat1<T>& x );

    int get_iterations() const { return its; }
};

#endif                          // __3T_MatVec_3T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/MatVec_3T.hh
//---------------------------------------------------------------------------//
