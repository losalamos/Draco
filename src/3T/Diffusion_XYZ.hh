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

using dsxx::Mat2;

#include "linalg/PCG_Ctrl.hh"
#include "linalg/Banded_Matrix.cc"

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
    Banded_Matrix< double, 7 > AA;

    PCG_Ctrl<double> pcg_ctrl;

    SP< MatVec_3T< MT, Diffusion_XYZ<MT> > > spmv;
    SP< PCG_PreCond<double> > precond;

    typename MT::cell_array vc;

    double dx, dy, dz;

  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::fcdsf FluxField;
    typedef typename MT::fcdsf DiscFluxField;
    
    typedef double NumT;

    Diffusion_XYZ( const SP<MT>& spm_, const pcg_DB& pcg_db );

    void solve( const typename MT::fcdsf& D,
		typename MT::cell_array& r,
		double dt, typename MT::cell_array& x,
		const typename MT::fcdsf& Eb );

    Mat2<double>& get_A() { return A; }
    SP< MatVec_3T< MT, Diffusion_XYZ<MT> > > get_matvec() { return spmv; }
    SP<MT> getMesh() { return spm; }
    const SP<MT> getMesh() const { return spm; }
};

#endif                          // __3T_Diffusion_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Diffusion_XYZ.hh
//---------------------------------------------------------------------------//
