//----------------------------------*-C++-*----------------------------------//
// Diffusion_XYZ.hh
// Geoffrey Furnish
// Thu Dec  4 09:30:51 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_Diffusion_XYZ_hh__
#define __3T_Diffusion_XYZ_hh__

#include "ds++/SP.hh"
#include "ds++/Mat.hh"

#include "linalg/PCG_Ctrl.hh"

#include "3T/MatVec_3T.hh"
#include "3T/PreCond.hh"

//===========================================================================//
// class Diffusion_XYZ - 

// 
//===========================================================================//

template<class MT>
class Diffusion_XYZ : private MT::Coord_Mapper {

    SP<MT> spm;

    Mat2<double> A;

    PCG_Ctrl<double> pcg_ctrl;

    SP< MatVec_3T< MT, Diffusion_XYZ<MT> > > spmv;
    SP< PCG_PreCond<double> > precond;

  public:
    typedef double NumT;

    Diffusion_XYZ( const SP<MT>& spm_, const pcg_DB& pcg_db );
    void solve( const typename MT::fcdsf& D,
		const typename MT::cell_array& rhs,
		double dt, typename MT::cell_array& x );

    Mat2<double>& get_A() { return A; }
};

#endif                          // __3T_Diffusion_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Diffusion_XYZ.hh
//---------------------------------------------------------------------------//
