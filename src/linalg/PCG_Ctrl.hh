//----------------------------------*-C++-*----------------------------------//
// PCG_Ctrl.hh
// Dave Nystrom
// Mon Jan 13 17:40:28 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __linalg_PCG_Ctrl_hh__
#define __linalg_PCG_Ctrl_hh__

#include "ds++/Mat.hh"
#include "ds++/SP.hh"

#include "linalg/pcg_DB.hh"
#include "linalg/PCG_MatVec.hh"
#include "linalg/PCG_PreCond.hh"

//===========================================================================//
// class PCG_Ctrl - 

// 
//===========================================================================//

template<class T>
class PCG_Ctrl : private pcg_DB
{
    Mat1<int>            iparm;
    Mat1<T>              fparm;
    Mat1<int>            iwk;
    Mat1<T>              fwk;
    Mat1<T>              xex;
    int                  nru;
    int                  ijob, ireq, iva, ivql, ivqr, ier;
    int                  itmeth, imatvec;

  public:
// Constructor.
    PCG_Ctrl( const pcg_DB& _pcg_db, int _nru );

// Public Methods
  public:
    void pcg_fe( Mat1<T>& x, const Mat1<T>& b,
		 SP< PCG_MatVec<T> > pcg_matvec,
		 SP< PCG_PreCond<T> > pcg_precond );

// Private Methods
  private:
    void set_default();
    void set_nwi();
    void set_nwf();
    void it_method( Mat1<T>& x, const Mat1<T>& b, Mat1<T>& xex );
    void print_params();
};

#endif                          // __linalg_PCG_Ctrl_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/PCG_Ctrl.hh
//---------------------------------------------------------------------------//
