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

//===========================================================================//
// class Diffusion_XYZ - 

// 
//===========================================================================//

template<class MT>
class Diffusion_XYZ : private MT::Coord_Mapper {

    SP<MT> spm;

    Mat2<double> A;

  public:
    Diffusion_XYZ( const SP<MT>& spm_ );
    void solve( const typename MT::fcdsf& D,
		const typename MT::cell_array& rhs,
		double dt, typename MT::cell_array& x );
};

#endif                          // __3T_Diffusion_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Diffusion_XYZ.hh
//---------------------------------------------------------------------------//
