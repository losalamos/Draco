//----------------------------------*-C++-*----------------------------------//
// PCG_PreCond.hh
// Dave Nystrom
// Fri May  9 12:30:54 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __linalg_PCG_PreCond_hh__
#define __linalg_PCG_PreCond_hh__

#include "ds++/Mat.hh"

//===========================================================================//
// class PCG_PreCond - 

// 
//===========================================================================//

template<class T>
class PCG_PreCond {
  public:
    virtual ~PCG_PreCond() {}

    virtual void  Left_PreCond( rtt_dsxx::Mat1<T>& x, const rtt_dsxx::Mat1<T>& b ) =0;
    virtual void Right_PreCond( rtt_dsxx::Mat1<T>& x, const rtt_dsxx::Mat1<T>& b ) =0;
};

#endif                          // __linalg_PCG_PreCond_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/PCG_PreCond.hh
//---------------------------------------------------------------------------//
