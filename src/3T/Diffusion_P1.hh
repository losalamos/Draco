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
    typedef typename MT::ccsf ccsf;
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::bssf bssf;

    dsxx::SP<MT> spm;

    Banded_Matrix< double, 7 > A;
    ccsf b;

    PCG_Ctrl<double> pcg_ctrl;

    dsxx::SP< MatVec_3T< Diffusion_P1<MT> > > spmv;
    dsxx::SP< PreCond< Diffusion_P1<MT> > > precond;

    dsxx::Mat1<double> dx, dy, dz;
    fcdsf fdeltal;

    fcdsf Dprime,      Dtwidle, Dhat;
    fcdsf Fprimeprime, Ftwidle, Fhat;

  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::fcdsf FluxField;
    typedef typename MT::fcdsf DiscFluxField;

    typedef double NumT;

    Diffusion_P1( const Diffusion_DB& diffdb,
                   const dsxx::SP<MT>& spm_, const pcg_DB& pcg_db );

    void solve( ccsf& phi,
                fcdsf& F,
		const fcdsf& D,
                const ccsf& sigmaabar,
                const ccsf& Qbar_r,
                const fcdsf& Fprime,
                const bssf& f_b );

    void solve( ccsf& phi,
		const fcdsf& D,
                const ccsf& sigmaabar,
                const ccsf& Qbar_r );

    void calculate_opposing_face( const fcdsf& X, fcdsf& Xprime );
    void calculate_Dtwidle( const fcdsf& D, const fcdsf& Dp );
    void calculate_Dhat_on_boundaries( const fcdsf& D );
    void calculate_Ftwidle( const fcdsf& D, const fcdsf& Fprime );
    void calculate_Fhat_on_boundaries( const fcdsf& D,
                                       const fcdsf& Fprime,
                                       const bssf& f_b );
    void calculate_A( const ccsf& sigmaabar );
    void calculate_b( const ccsf& Qbar_r );
    void calculate_b( const ccsf& Qbar_r, const fcdsf& Fh );

    void solve_A_phi_equals_b( ccsf& phi );

    void calculate_new_F( const fcdsf& D, const fcdsf& Fprime,
                          const bssf& f_b, const ccsf& phi,
			  fcdsf& F );

    Banded_Matrix<double,7>& get_A() { return A; }
};

#endif                          // __3T_Diffusion_P1_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/Diffusion_P1.hh
//---------------------------------------------------------------------------//
