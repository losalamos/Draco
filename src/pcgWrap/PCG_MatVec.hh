//----------------------------------*-C++-*----------------------------------//
// PCG_MatVec.hh
// Dave Nystrom
// Fri May  9 12:30:20 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __pcgWrap_PCG_MatVec_hh__
#define __pcgWrap_PCG_MatVec_hh__

#include "ds++/Mat.hh"

//===========================================================================//
// class PCG_MatVec - 

// 
//===========================================================================//

namespace rtt_pcgWrap {

template<class T>
class PCG_MatVec {
  public:
    virtual ~PCG_MatVec() {}

    virtual void MatVec( rtt_dsxx::Mat1<T>& b,
			 const rtt_dsxx::Mat1<T>& x ) = 0;
};

} // namespace rtt_pcgWrap

#endif                          // __pcgWrap_PCG_MatVec_hh__

//---------------------------------------------------------------------------//
//                              end of pcgWrap/PCG_MatVec.hh
//---------------------------------------------------------------------------//
