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

#include "linalg/PCG_Ctrl.hh"
#include "linalg/Banded_Matrix.cc"

#include "3T/Diffusion_DB.hh"
#include "3T/MatVec_3T.hh"
#include "3T/PreCond.hh"

//===========================================================================//
// class Diffusion_XYZ - 

// 
//===========================================================================//

template<class MT>
class Diffusion_XYZ : private MT::Coord_Mapper,
                      protected Diffusion_DB
{
    SP<MT> spm;

    Banded_Matrix< double, 7 > A;

    PCG_Ctrl<double> pcg_ctrl;

    SP< MatVec_3T< Diffusion_XYZ<MT> > > spmv;
    SP< PreCond< Diffusion_XYZ<MT> > > precond;
 
    typedef typename MT::ccsf cell_array_double;
//     typename MT::template cell_array<double> vc;
    cell_array_double vc;

    double dx, dy, dz;

  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::fcdsf FluxField;
    typedef typename MT::fcdsf DiscFluxField;
    
    typedef double NumT;

    Diffusion_XYZ( const Diffusion_DB& diffdb,
                   const SP<MT>& spm_, const pcg_DB& pcg_db );

    void solve( const typename MT::fcdsf& D,
		cell_array_double& r,
		double dt,
                cell_array_double& x,
		const typename MT::fcdsf& Eb );

    Banded_Matrix<double,7>& get_A() { return A; }
    SP< MatVec_3T< Diffusion_XYZ<MT> > > get_matvec() { return spmv; }
    SP<MT> getMesh() { return spm; }
    const SP<MT> getMesh() const { return spm; }
};

#endif                          // __3T_Diffusion_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Diffusion_XYZ.hh
//---------------------------------------------------------------------------//
