//----------------------------------*-C++-*----------------------------------//
// dummy_package.cc
// John McGhee
// Thu Aug 27 07:48:41 1998
//---------------------------------------------------------------------------//
// @> A dummy package to exercize the time-step controller field advisors.
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>
#include "ds++/Soft_Equivalence.hh"
#include "c4/global.hh"
#include "../ts_manager.hh"
#include "../field_ts_advisor.hh"
#include "dummy_package.hh"

using namespace rtt_timestep;

std::vector<double> operator*(double lhs, const std::vector<double> &rhs)
{
    std::vector<double> results(rhs.size());

    for (int i=0; i<rhs.size(); i++)
	results[i] = lhs*rhs[i];

    return results;
}

dummy_package::dummy_package(ts_manager &tsm_)
    : tsm(tsm_)
{

// Set up a Electron Temperature advisor

    sp_te = new field_ts_advisor( "Electron Temperature",
				   ts_advisor::max,
				   field_ts_advisor::a_mean);
    tsm.add_advisor(sp_te);

// Set up a Ion-Temperature advisor

    sp_ti = new field_ts_advisor("Ion Temperature",
				  ts_advisor::max,
				  field_ts_advisor::rc_mean);
    tsm.add_advisor(sp_ti);

// Set up a Radiation-Intensity advisor

    sp_ri = new field_ts_advisor("Radiation Intensity",
				  ts_advisor::max,
				  field_ts_advisor::rcq_mean);
    tsm.add_advisor(sp_ri);
}

dummy_package::~dummy_package()
{
    tsm.remove_advisor(sp_te);
    tsm.remove_advisor(sp_ti);
    tsm.remove_advisor(sp_ri);
}

void dummy_package::advance_state()
{
    // RMR The following code does not adequately cover
    // all of the possible failure modes.
    // Please update this method with fuller coverage.

    // Create a set of dummy arrays to serve as control fields for
    // use in exercizing the various advisors.

    const double a1Array[] = {1., 10., 11., 3., 2., 5., 5., 6.7};
    const int sizeaArray = sizeof(a1Array)/sizeof(a1Array[0]);

    const int sizea = sizeaArray;
    const double *a1 = a1Array;
    
    std::vector<double> te_old(a1,a1+sizea);
    std::vector<double> te_new = 1.09*te_old;
    std::vector<double> ti_old=0.97*te_old;
    std::vector<double> ti_new=1.05*te_old;
    std::vector<double> ri_old=1.10*te_old;
    std::vector<double> ri_new=1.15*te_old;

// Set a floor for the electron temperature controller, to
// execcize this method. Just accelpt the default floor on the
// other controllers.

    sp_te -> set_floor(te_new,0.001);

// Get a new time-step from each of the advisors that
// belong to this package.

    sp_te -> update_tstep(tsm, te_old, te_new);
    sp_ti -> update_tstep(tsm, ti_old, ti_new);
    sp_ri -> update_tstep(tsm, ri_old, ri_new);

    return;
}

bool dummy_package::tests_passed() const
{
    using rtt_dsxx::soft_equiv;

    double const prec( 1.0e-5 );
    double const ref1( 1.371742 );
    double const ref2( 1.496914 );
    double const ref3( 2.716049 );
    
    return soft_equiv( ref1, sp_te->get_dt_rec(tsm), prec )
	&& soft_equiv( ref2, sp_ti->get_dt_rec(tsm), prec )
	&& soft_equiv( ref3, sp_ri->get_dt_rec(tsm), prec );
}

//---------------------------------------------------------------------------//
//                              end of dummy_package.cc
//---------------------------------------------------------------------------//
