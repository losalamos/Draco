//----------------------------------*-C++-*----------------------------------//
// PCG_PreCond.hh
// Dave Nystrom
// Fri May  9 12:30:54 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __pcgWrap_PCG_PreCond_hh__
#define __pcgWrap_PCG_PreCond_hh__

#include "ds++/Mat.hh"

//===========================================================================//
// class PCG_PreCond - 

// 
//===========================================================================//

namespace rtt_pcgWrap {

template<class T>
class PCG_PreCond {
  public:
    virtual ~PCG_PreCond() {}

    virtual void Left_PreCond( rtt_dsxx::Mat1<T>& x,
			       const rtt_dsxx::Mat1<T>& b ) = 0;
    
    virtual void Right_PreCond( rtt_dsxx::Mat1<T>& x,
				const rtt_dsxx::Mat1<T>& b ) = 0;
};

} // namespace rtt_pcgWrap

#endif                          // __pcgWrap_PCG_PreCond_hh__

//---------------------------------------------------------------------------//
//                              end of pcgWrap/PCG_PreCond.hh
//---------------------------------------------------------------------------//
