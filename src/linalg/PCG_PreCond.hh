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

using dsxx::Mat1;

//===========================================================================//
// class PCG_PreCond - 

// 
//===========================================================================//

template<class T>
class PCG_PreCond {

  public:
//     PCG_PreCond();
    virtual ~PCG_PreCond() {}

    virtual void  Left_PreCond( Mat1<T>& x, Mat1<T>& b ) =0;
    virtual void Right_PreCond( Mat1<T>& x, Mat1<T>& b ) =0;
};

#endif                          // __linalg_PCG_PreCond_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/PCG_PreCond.hh
//---------------------------------------------------------------------------//
