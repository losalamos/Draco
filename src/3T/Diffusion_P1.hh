//----------------------------------*-C++-*----------------------------------//
// Diffusion_P1.hh
// Geoffrey M. Furnish
// Thu May 28 13:16:55 1998
//---------------------------------------------------------------------------//
// @> Class for solving P1 Diffusion equations to support P1-3T package.
//---------------------------------------------------------------------------//

#ifndef __3T_Diffusion_P1_hh__
#define __3T_Diffusion_P1_hh__

#include "ds++/SP.hh"

#include "linalg/PCG_Ctrl.hh"
#include "linalg/Banded_Matrix.hh"

#include "3T/Diffusion_DB.hh"
#include "3T/MatVec_3T.hh"
#include "3T/PreCond.hh"

//===========================================================================//
// class Diffusion_P1 - Solve P1 diffusion equation for P13T package

// This class implements the P1 diffusion solver needed to support the
// conduction and radiation diffusion equations in the P13T package.  These
// equations are documented in "3-T Diffusion with Material Motion
// Corrections", part of the Draco doc set.
//===========================================================================//

template<class MT>
class Diffusion_P1 : private MT::Coord_Mapper,
                     protected Diffusion_DB
{
    typedef typename MT::template cell_array<double> ccsf;
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::template bssf<double> bssf;

    SP<MT> spm;

    Banded_Matrix< double, 7 > A;

    PCG_Ctrl<double> pcg_ctrl;

    SP< MatVec_3T< Diffusion_P1<MT> > > spmv;
    SP< PreCond< Diffusion_P1<MT> > > precond;

    Mat1<double> dx, dy, dz;

    fcdsf Dprime, Dtwidle, Dhat;

  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::fcdsf FluxField;
    typedef typename MT::fcdsf DiscFluxField;

    typedef double NumT;

    Diffusion_P1( const Diffusion_DB& diffdb,
                   const SP<MT>& spm_, const pcg_DB& pcg_db );

    void solve( const fcdsf& D,
                const ccsf& sigmaabar,
                const ccsf& Qbar_r,
                const fcdsf& Fprime,
                const bssf& f_b,
                ccsf& phi,
                fcdsf& F );

    void solve( const fcdsf& D,
                const ccsf& sigmaabar,
                const ccsf& Qbar_r,
                ccsf& phi );

    void calculate_Dprime( const fcdsf& D );
    void calculate_Dtwidle( const fcdsf& D, const fcdsf& Dp );
    void calculate_Dhat_on_boundaries();

    Banded_Matrix<double,7>& get_A() { return A; }
};

#endif                          // __3T_Diffusion_P1_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Diffusion_P1.hh
//---------------------------------------------------------------------------//
