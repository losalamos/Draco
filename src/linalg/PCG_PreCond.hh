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
    PCG_PreCond();
    ~PCG_PreCond();

    virtual void  Left_PreCond( Mat1<T>& x, Mat1<T>& b );
    virtual void Right_PreCond( Mat1<T>& x, Mat1<T>& b );
};

#endif                          // __linalg_PCG_PreCond_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/PCG_PreCond.hh
//---------------------------------------------------------------------------//
