//----------------------------------*-C++-*----------------------------------//
// TstPCG_MatVec.hh
// Dave Nystrom
// Fri May  9 13:39:23 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __linalg_test_TstPCG_MatVec_hh__
#define __linalg_test_TstPCG_MatVec_hh__

#include "linalg/PCG_MatVec.hh"

#include "ds++/Mat.hh"

//===========================================================================//
// class TstPCG_MatVec - 

// 
//===========================================================================//

template<class T>
class TstPCG_MatVec : public PCG_MatVec<T> {
    int nxs, nys;

  public:
    TstPCG_MatVec( int _nxs, int _nys );
    ~TstPCG_MatVec();

    void MatVec( Mat1<T>& b, Mat1<T>& x );
};

#endif                          // __linalg_test_TstPCG_MatVec_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/test/TstPCG_MatVec.hh
//---------------------------------------------------------------------------//
