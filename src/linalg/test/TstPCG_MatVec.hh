//----------------------------------*-C++-*----------------------------------//
// TstPCG_MatVec.hh
// Dave Nystrom
// Fri May  9 13:39:23 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __pcgWrap_test_TstPCG_MatVec_hh__
#define __pcgWrap_test_TstPCG_MatVec_hh__

#include "../PCG_MatVec.hh"

#include "ds++/Mat.hh"

//===========================================================================//
// class TstPCG_MatVec - 

// 
//===========================================================================//

template<class T>
class TstPCG_MatVec : public rtt_pcgWrap::PCG_MatVec<T>
{
    int nxs, nys;

  public:
    TstPCG_MatVec( int _nxs, int _nys );
    ~TstPCG_MatVec();

    void MatVec( rtt_dsxx::Mat1<T>& b,
		 const rtt_dsxx::Mat1<T>& x );
};

#endif                          // __pcgWrap_test_TstPCG_MatVec_hh__

//---------------------------------------------------------------------------//
//                              end of pcgWrap/test/TstPCG_MatVec.hh
//---------------------------------------------------------------------------//
