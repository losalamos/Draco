//----------------------------------*-C++-*----------------------------------//
// Inner_iter.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Controller of inner iterations.
//---------------------------------------------------------------------------//

#include "sn/Inner_iter.hh"
#include "sn/timer.hh"

void Inner_iter::do_inner_iter( Input_edit &data )
{
    // This routine performs the inner iteration sweeps.
    // It also tests for convergence, does timing and prints results.

    // initialize all flux moments and the iteration counter to zero
      
    Mat4<REAL> flux(data.jt(),data.kt(),data.it(),data.nm());

    int its = 0;

    // begin cpu and wallclock timing

    REAL cpu0;   // starting time for determining total cpu       time
    REAL wall0;  // starting time for determining total wallclock time

    timer( cpu0, wall0 );

    // instantiate a Cross_section object

    Cross_section xsec( data );

    // initialize source

    Source src( data );

    Mat4<REAL> src_mom(data.jt(),data.kt(),data.it(),data.nm());  // src moments
                                                                  // per cell 
    // initialize error

    Error converge( data ) ;

    // instantiate an Sn_constants object

    Sn_constants sn( data );

    // pre-calculate all possible factors for use in inner solver loops

    Pre_calcs pre( data, xsec, sn );

    // begin iterations

    Mat1<REAL> lkgs_l(6);  // leakages from each face of the cube

    while (1)
    {

        its += 1;

        // compute the full source, all moments

        src.build_source( data, xsec, src_mom, flux );

        // initialize the error convergence check

        converge.error_init( data, flux );

        // re-initialize the current flux moments to zero

        flux = 0.0;

        // perform an ordered sweep through the mesh.

        Sweep3d sweep3d_object;

        sweep3d_object.do_sweep( data, xsec, sn, pre, lkgs_l, src_mom, flux );

        // compute the convergence error

        converge.error_calc( data, flux );
      
        // print iteration/error information and test for convergence

        converge.error_print( its );

        if ( converge.error_test( data, its ) )
            break;

    }

    // stop cpu and wallclock timing

    REAL cpu1;   // ending   time for determining total cpu       time
    REAL wall1;  // ending   time for determining total wallclock time

    timer( cpu1, wall1 );

    // perform edits and outputs

    Output_edit output;

    output.print_flux_sum( data, flux );
    output.print_leakages( lkgs_l );
    output.print_timing( cpu0, cpu1, wall0, wall1, data, its );
    output.print_flux_moments( data, flux );

    return;
}

//---------------------------------------------------------------------------//
//                              end of Inner_iter.cc
//---------------------------------------------------------------------------//

